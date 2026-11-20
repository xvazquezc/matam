#include "perfectOverlap.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <stdint.h>
#include <vector>

// The vectors below only live inside this file; GCC warns that returning
// them would change the ABI on targets without AVX.
#pragma GCC diagnostic ignored "-Wpsabi"

namespace
{

// 16 cells of an anti-diagonal at a time. GCC vector extensions compile to
// whatever SIMD the target has (two SSE2 registers on plain x86-64).
int const LANES = 16;
typedef int16_t CellVector __attribute__((vector_size(LANES * sizeof(int16_t))));
typedef uint8_t BaseVector __attribute__((vector_size(LANES)));
typedef int8_t BaseMaskVector __attribute__((vector_size(LANES)));

// Scores stay within +-(length I + length J), lanes past the end included.
std::size_t const MAX_LENGTH_SUM = std::numeric_limits<int16_t>::max() - 2 * LANES;

struct Workspace
{
    std::vector<uint8_t> basesI;         // read i
    std::vector<uint8_t> reversedBasesJ; // read j, reversed
    std::vector<int16_t> diagonals[3];   // anti-diagonals k, k-1, k-2 by row i
};

thread_local Workspace workspace;

inline CellVector loadCells(int16_t const *cells)
{
    CellVector v;
    std::memcpy(&v, cells, sizeof v);
    return v;
}

inline BaseVector loadBases(uint8_t const *bases)
{
    BaseVector v;
    std::memcpy(&v, bases, sizeof v);
    return v;
}

inline CellVector maxCells(CellVector a, CellVector b)
{
    return a > b ? a : b;
}

/******************************************************************************
    Best score of the overlap alignment (+1 match, -1 mismatch or gap, free
    end gaps) and where it is reached on the last row or column.

    H[i][j] = max(H[i-1][j-1] + s(i, j), H[i-1][j] - 1, H[i][j-1] - 1), with
    H[0][j] = H[i][0] = 0. Cells of an anti-diagonal k = i + j only depend on
    anti-diagonals k-1 and k-2, so they are computed LANES at a time.
******************************************************************************/
struct BestEnd
{
    int score;
    int count;  // number of cells of the last row and column reaching score
    int i;
    int j;
};

BestEnd bestOverlapEnd(seqan::String<seqan::Iupac> const &readSeqI,
                       seqan::String<seqan::Iupac> const &readSeqJ)
{
    int const n = seqan::length(readSeqI);
    int const m = seqan::length(readSeqJ);

    Workspace &w = workspace;
    w.basesI.assign(n + LANES, 0);
    w.reversedBasesJ.assign(m + LANES, 0);
    for (int t = 0; t < n; ++t)
        w.basesI[t] = seqan::ordValue(readSeqI[t]);
    for (int t = 0; t < m; ++t)
        w.reversedBasesJ[t] = seqan::ordValue(readSeqJ[m - 1 - t]);
    for (auto &diagonal : w.diagonals)
        diagonal.assign(n + 1 + LANES, 0);

    // Anti-diagonals 0 and 1 hold boundary cells only, which are all 0.
    int16_t *previous2 = w.diagonals[0].data();
    int16_t *previous1 = w.diagonals[1].data();
    int16_t *current = w.diagonals[2].data();

    BestEnd best = {std::numeric_limits<int>::min(), 0, 0, 0};
    auto consider = [&best](int score, int i, int j) {
        if (score > best.score)
            best = {score, 1, i, j};
        else if (score == best.score)
            ++best.count;
    };

    uint8_t const *basesI = w.basesI.data();
    uint8_t const *reversedBasesJ = w.reversedBasesJ.data();
    CellVector const ones = CellVector{} + 1;

    for (int k = 2; k <= n + m; ++k)
    {
        std::swap(previous2, previous1);  // previous2 <- k-2, previous1 <- old k-3
        std::swap(previous1, current);    // previous1 <- k-1, current <- old k-3

        int const lo = std::max(1, k - m);
        int const hi = std::min(n, k - 1);
        // Cell (i, j = k - i) compares base i-1 of read i with base j-1 of
        // read j, which is reversedBasesJ[m - k + i].
        uint8_t const *rowBasesJ = reversedBasesJ + (m - k);
        for (int i = lo; i <= hi; i += LANES)
        {
            CellVector up = loadCells(previous1 + i - 1);   // (i-1, j)
            CellVector left = loadCells(previous1 + i);     // (i, j-1)
            CellVector diag = loadCells(previous2 + i - 1); // (i-1, j-1)
            BaseMaskVector equal = loadBases(basesI + i - 1) == loadBases(rowBasesJ + i);
            // equal lanes are -1, others 0: substitution score +1 or -1
            CellVector substitution = -ones - 2 * __builtin_convertvector(equal, CellVector);
            CellVector cells = maxCells(diag + substitution, maxCells(up, left) - ones);
            std::memcpy(current + i, &cells, sizeof cells);
        }
        // Boundary cells; set after the loop, whose last lanes may run past hi.
        if (k <= m)
            current[0] = 0;
        if (k <= n)
            current[k] = 0;

        if (k - n >= 1 && k - n <= m)          // last row: (n, k - n)
            consider(current[n], n, k - n);
        if (k - m >= 1 && k - m < n)           // last column: (k - m, m)
            consider(current[k - m], k - m, m);
    }
    return best;
}

/******************************************************************************
    Source position of a view position in an ungapped row placed at
    viewBegin, as SeqAn's toSourcePosition gives it: 0 in the leading gaps,
    the sequence length in the trailing gaps.
******************************************************************************/
inline int sourcePosition(int viewPos, int viewBegin, int length)
{
    if (viewPos < viewBegin)
        return 0;
    return std::min(viewPos - viewBegin, length);
}

}  // namespace

bool computePerfectOverlap(OverlapStatistics &alignOvStats,
                           OverlapDescription &alignOvDescription,
                           seqan::String<seqan::Iupac> const &readSeqI,
                           seqan::String<seqan::Iupac> const &readSeqJ)
{
    static_assert(sizeof(seqan::Iupac) == 1, "one byte per base");

    std::size_t const n = seqan::length(readSeqI);
    std::size_t const m = seqan::length(readSeqJ);
    if (n == 0 || m == 0 || n + m > MAX_LENGTH_SUM)
        return false;

    BestEnd const best = bestOverlapEnd(readSeqI, readSeqJ);
    int const diagonalLength = std::min(best.i, best.j);
    if (best.count != 1 || best.score != diagonalLength)
        return false;

    // The alignment is the ungapped diagonal ending at cell (best.i, best.j):
    // base x of read i faces base x - offset of read j. Laid out as SeqAn
    // does, each read starts after the leading gaps of the other.
    int const offset = best.i - best.j;
    int const lengthI = n;
    int const lengthJ = m;
    int const lastPosI = lengthI - 1;
    int const lastPosJ = lengthJ - 1;
    int const beginPosInAlignI = std::max(0, -offset);
    int const beginPosInAlignJ = std::max(0, offset);
    int const endPosInAlignI = beginPosInAlignI + lastPosI;
    int const endPosInAlignJ = beginPosInAlignJ + lastPosJ;

    // From here on, the same computations as computeAlignmentStats.
    alignOvDescription.readLengthI = lengthI;
    alignOvDescription.readLengthJ = lengthJ;

    int const overlapBeginPos = std::max(beginPosInAlignI, beginPosInAlignJ);
    int const overlapEndPos = std::min(endPosInAlignI, endPosInAlignJ);
    alignOvStats.totalOverlapPositions = overlapEndPos - overlapBeginPos + 1;

    alignOvDescription.ovBeginPosI = std::max(0, sourcePosition(beginPosInAlignJ, beginPosInAlignI, lengthI));
    alignOvDescription.ovBeginPosJ = std::max(0, sourcePosition(beginPosInAlignI, beginPosInAlignJ, lengthJ));

    alignOvDescription.ovEndPosI = std::min(lastPosI, sourcePosition(endPosInAlignJ, beginPosInAlignI, lengthI));
    alignOvDescription.ovEndPosJ = std::min(lastPosJ, sourcePosition(endPosInAlignI, beginPosInAlignJ, lengthJ));

    if ((alignOvDescription.ovBeginPosI == 0 && alignOvDescription.ovEndPosI == lastPosI)
        || (alignOvDescription.ovBeginPosJ == 0 && alignOvDescription.ovEndPosJ == lastPosJ))
    {
        alignOvDescription.isEnclosedOverlap = true;
    }

    // Columns of an enclosed overlap are not counted, as in computeAlignmentStats.
    if (!alignOvDescription.isEnclosedOverlap)
        alignOvStats.matchesNum = alignOvStats.totalOverlapPositions;

    return true;
}

bool sameOverlap(OverlapStatistics const &statsA,
                 OverlapDescription const &descriptionA,
                 OverlapStatistics const &statsB,
                 OverlapDescription const &descriptionB)
{
    return statsA.totalOverlapPositions == statsB.totalOverlapPositions
        && statsA.coherentPositions == statsB.coherentPositions
        && statsA.indelNum == statsB.indelNum
        && statsA.matchesNum == statsB.matchesNum
        && statsA.mismatchesNum == statsB.mismatchesNum
        && descriptionA.readIDI == descriptionB.readIDI
        && descriptionA.readIDJ == descriptionB.readIDJ
        && descriptionA.readLengthI == descriptionB.readLengthI
        && descriptionA.readLengthJ == descriptionB.readLengthJ
        && descriptionA.ovBeginPosI == descriptionB.ovBeginPosI
        && descriptionA.ovEndPosI == descriptionB.ovEndPosI
        && descriptionA.ovBeginPosJ == descriptionB.ovBeginPosJ
        && descriptionA.ovEndPosJ == descriptionB.ovEndPosJ
        && descriptionA.isEnclosedOverlap == descriptionB.isEnclosedOverlap
        && descriptionA.isRead2RC == descriptionB.isRead2RC;
}
