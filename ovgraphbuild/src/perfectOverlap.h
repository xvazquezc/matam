#ifndef PERFECTOVERLAP_H
#define PERFECTOVERLAP_H

#include <seqan/sequence.h>
#include <seqan/modifier.h>

#include "alignmentsComparison.h"
#include "assemblyStatistics.h"

/******************************************************************************
    Shortcut for the de novo overlap alignment of two reads

    The re-alignment of two compatible reads (align2ReadSequences followed by
    computeAlignmentStats) costs a full SeqAn alignment with traceback. Most of
    these reads overlap without any difference, and for them the alignment
    SeqAn returns can be proven to be that ungapped overlap:

    - The overlap alignment scores +1 per match and -1 per mismatch or gap,
      with free end gaps. Take the cell E of the last row or column where an
      alignment ends, and L the length of the diagonal from E back to the
      first row or column. Every alignment ending at E scores at most L, and
      only the ungapped all-match diagonal scores L: any other path has fewer
      aligned pairs or pays for its gaps.
    - So when the best score over the last row and column is reached at one
      cell only, and equals L there, SeqAn has a single optimal alignment to
      return: that diagonal. No tie-breaking is involved.

    computePerfectOverlap computes the alignment score only (vectorized over
    the anti-diagonals of the matrix) and checks these two conditions. When
    they hold, it fills the statistics and overlap description exactly as
    computeAlignmentStats would and returns true. Otherwise it returns false
    and the caller runs the SeqAn alignment.
******************************************************************************/
bool computePerfectOverlap(OverlapStatistics &,
                           OverlapDescription &,
                           seqan::String<seqan::Iupac> const &, // read i
                           seqan::String<seqan::Iupac> const &); // read j

/******************************************************************************
    Whether two re-alignments gave the same statistics and overlap
******************************************************************************/
bool sameOverlap(OverlapStatistics const &,
                 OverlapDescription const &,
                 OverlapStatistics const &,
                 OverlapDescription const &);

#endif // PERFECTOVERLAP_H
