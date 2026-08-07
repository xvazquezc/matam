#include "compatibilityGraphBuilding.h"

#include <atomic>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <thread>

namespace
{
constexpr int MAX_WORKER_THREADS = 24;
constexpr int CHUNKS_PER_WORKER = 8;

struct PairRange
{
    int64_t begin;
    int64_t end;
    int64_t pairCount;
};

struct PairChunkResult
{
    PairRange range;
    GlobalStatistics statistics;
    int64_t completedPairs = 0;
    std::string asqgFragment;
    std::string csvEdgesFragment;
    std::string error;
};

int64_t totalPairs(int64_t readCount)
{
    return readCount < 2 ? 0 : readCount * (readCount - 1) / 2;
}

int64_t pairsBefore(int64_t readCount, int64_t outerIndex)
{
    return outerIndex * (2 * readCount - outerIndex - 1) / 2;
}

int64_t balancedTarget(int64_t begin,
                       int64_t end,
                       int partIndex,
                       int partCount)
{
    int64_t size = end - begin;
    int64_t quotient = size / partCount;
    int64_t remainder = size % partCount;
    return begin + quotient * partIndex + std::min<int64_t>(partIndex, remainder);
}

int64_t nearestOuterBoundary(int64_t readCount,
                             int64_t minOuterIndex,
                             int64_t maxOuterIndex,
                             int64_t targetPairCount)
{
    int64_t low = minOuterIndex;
    int64_t high = maxOuterIndex;

    while (low < high)
    {
        int64_t mid = low + (high - low) / 2;
        if (pairsBefore(readCount, mid) < targetPairCount)
            low = mid + 1;
        else
            high = mid;
    }

    if (low > minOuterIndex)
    {
        int64_t previous = low - 1;
        int64_t previousDistance = targetPairCount - pairsBefore(readCount, previous);
        int64_t currentDistance = pairsBefore(readCount, low) - targetPairCount;
        if (previousDistance <= currentDistance)
            return previous;
    }

    return low;
}

PairRange partitionPairRange(int64_t readCount,
                             int64_t outerBegin,
                             int64_t outerEnd,
                             int partIndex,
                             int partCount)
{
    int64_t pairBegin = pairsBefore(readCount, outerBegin);
    int64_t pairEnd = pairsBefore(readCount, outerEnd);
    int64_t beginTarget = balancedTarget(pairBegin, pairEnd, partIndex, partCount);
    int64_t endTarget = balancedTarget(pairBegin, pairEnd, partIndex + 1, partCount);

    int64_t begin = nearestOuterBoundary(readCount, outerBegin, outerEnd, beginTarget);
    int64_t end = nearestOuterBoundary(readCount, begin, outerEnd, endTarget);

    return PairRange{begin, end, pairsBefore(readCount, end) - pairsBefore(readCount, begin)};
}

void mergeCompatibilityStatistics(GlobalStatistics &destination,
                                  GlobalStatistics const &source)
{
    destination.numConsideredReadsPairs += source.numConsideredReadsPairs;
    destination.readPairsWithCommonRefNum += source.readPairsWithCommonRefNum;
    destination.numConsideredAlignmentsPairs += source.numConsideredAlignmentsPairs;
    destination.alignmentsPairsWithCommonRefNum += source.alignmentsPairsWithCommonRefNum;
    destination.numPotentialOverlaps += source.numPotentialOverlaps;
    destination.numOverlaps += source.numOverlaps;
    destination.enclosedOverlapsNum += source.enclosedOverlapsNum;
    destination.compatibleReadPairsNum += source.compatibleReadPairsNum;
    destination.incompatibleReadPairsNum += source.incompatibleReadPairsNum;
    destination.neitherCompatNorIncompatReadPairsNum += source.neitherCompatNorIncompatReadPairsNum;
    destination.truePositive += source.truePositive;
    destination.falsePositive += source.falsePositive;
    destination.trueNegative += source.trueNegative;
    destination.falseNegative += source.falseNegative;
}

bool appendFile(std::ostream &destination, std::string const &sourceFilename)
{
    std::ifstream source(sourceFilename, std::ifstream::in);
    if (!source)
        return false;

    source.seekg(0, std::ios::end);
    if (source.tellg() == 0)
        return true;
    source.seekg(0, std::ios::beg);

    destination << source.rdbuf();
    return destination.good();
}

void removeFragments(std::vector<PairChunkResult> const &chunkResults)
{
    for (auto const &result : chunkResults)
    {
        if (!result.asqgFragment.empty())
            std::remove(result.asqgFragment.c_str());
        if (!result.csvEdgesFragment.empty())
            std::remove(result.csvEdgesFragment.c_str());
    }
}
}

//template<typename T> class TD;

/******************************************************************************
    Build the compatibility graph
******************************************************************************/
void buildCompatibilityGraph(TGraph &graph,
                             std::vector<TVertexDescriptor > &vertices,
                             TProperties &readNames,
                             GlobalStatistics &globalStats,
                             std::vector<std::vector<seqan::BamAlignmentRecord> > const &bamRecordBuffer,
                             AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    // Initialize graph vertices and store read names
    initializeGraph(graph, vertices, readNames, bamRecordBuffer, options);

    // Compute the compatibility for the read pairs assigned to this shard.
    int64_t mappedReadsNum = bamRecordBuffer.size();
    int64_t maxReadsPairs = totalPairs(mappedReadsNum);
    int64_t outerEnd = std::max<int64_t>(0, mappedReadsNum - 1);
    PairRange shardRange = partitionPairRange(mappedReadsNum,
                                              0,
                                              outerEnd,
                                              options.pairShardIndex,
                                              options.pairShardCount);

    if (options.debug)
    {
        std::cout << "DEBUG: Computing compatibility graph"
                  << "\n" << "\n" << std::flush;
    }

    // for performance reason to improve stream speed
    std::ios::sync_with_stdio(false);

    // Declare output files for the overlap graph
    std::ofstream asqgFile;
    std::ofstream csvNodesFile;
    std::ofstream csvEdgesFile;

    // Initialise the ASQG output file if needed
    if (options.outputASQG)
    {
        std::string asqgFilename(seqan::toCString(options.outputBasename));
        asqgFilename += ".asqg";

        // Opening ASQG output file
        asqgFile.open(asqgFilename, std::ofstream::out | std::ofstream::trunc);
        if (!asqgFile)
            throw std::runtime_error("could not open ASQG output file");

        // Write ASQG header
        asqgFile << "HT\tVN:i:1\tER:f:0\t"
                 << "OL:i:" << options.minOverlapLength << "\t"
                 << "IN:Z:" << options.mySamFile << "\t"
                 << "CN:i:0\tTE:i:0\n";

        // Write ASQG nodes (aka reads)
        writeReadsToASQG(asqgFile, bamRecordBuffer, options);
    }

    // Initialise the CSV output files if needed
    if (options.outputCSV)
    {
        std::string csvNodesFilename(seqan::toCString(options.outputBasename));
        csvNodesFilename += ".nodes.csv";
        std::string csvEdgesFilename(seqan::toCString(options.outputBasename));
        csvEdgesFilename += ".edges.csv";

        // Opening CSV output files
        csvNodesFile.open(csvNodesFilename, std::ofstream::out | std::ofstream::trunc);
        csvEdgesFile.open(csvEdgesFilename, std::ofstream::out | std::ofstream::trunc);
        if (!csvNodesFile || !csvEdgesFile)
            throw std::runtime_error("could not open CSV output files");

        // Write CSV nodes
        writeReadsToCSV(csvNodesFile, readNames, bamRecordBuffer, options);
        csvNodesFile.close();

        // Write CSV Edges header
        csvEdgesFile << "Source;Target;Type;Coemitted;Weight;Pid;MultiRef" << "\n";
    }

    int cappedThreadCount = std::min(options.threads, MAX_WORKER_THREADS);
    if (options.threads > MAX_WORKER_THREADS)
    {
        std::cout << "INFO: Capping graph worker threads at "
                  << MAX_WORKER_THREADS << " (requested " << options.threads << ")\n";
    }

    int64_t assignedOuterIndices = shardRange.end - shardRange.begin;
    int targetChunkCount = cappedThreadCount == 1
        ? 1
        : cappedThreadCount * CHUNKS_PER_WORKER;
    targetChunkCount = std::min<int64_t>(targetChunkCount,
                                         std::max<int64_t>(1, assignedOuterIndices));

    std::vector<PairChunkResult> chunkResults;
    chunkResults.reserve(targetChunkCount);
    for (int chunkIndex = 0; chunkIndex < targetChunkCount; ++chunkIndex)
    {
        PairRange range = partitionPairRange(mappedReadsNum,
                                             shardRange.begin,
                                             shardRange.end,
                                             chunkIndex,
                                             targetChunkCount);
        if (range.pairCount == 0)
            continue;

        chunkResults.emplace_back();
        chunkResults.back().range = range;
    }

    int workerCount = std::min<int>(cappedThreadCount, chunkResults.size());
    std::vector<std::thread> workers;
    std::string outputBasename(seqan::toCString(options.outputBasename));
    std::atomic<size_t> nextChunk(0);
    std::atomic<bool> stopWorkers(false);

    try
    {
        for (int workerIndex = 0; workerIndex < workerCount; ++workerIndex)
        {
            workers.emplace_back([&]()
            {
                while (!stopWorkers.load())
                {
                    size_t chunkIndex = nextChunk.fetch_add(1);
                    if (chunkIndex >= chunkResults.size())
                        break;

                    PairChunkResult &result = chunkResults[chunkIndex];
                    std::ofstream chunkAsqgFile;
                    std::ofstream chunkCsvEdgesFile;

                    try
                    {
                        if (options.outputASQG)
                        {
                            result.asqgFragment = outputBasename + ".asqg.chunk-"
                                                  + std::to_string(chunkIndex) + ".tmp";
                            chunkAsqgFile.open(result.asqgFragment,
                                               std::ofstream::out | std::ofstream::trunc);
                            if (!chunkAsqgFile)
                                throw std::runtime_error("could not open ASQG chunk fragment");
                        }

                        if (options.outputCSV)
                        {
                            result.csvEdgesFragment = outputBasename + ".edges.csv.chunk-"
                                                      + std::to_string(chunkIndex) + ".tmp";
                            chunkCsvEdgesFile.open(result.csvEdgesFragment,
                                                   std::ofstream::out | std::ofstream::trunc);
                            if (!chunkCsvEdgesFile)
                                throw std::runtime_error("could not open CSV chunk fragment");
                        }

                        for (int64_t i = result.range.begin; i < result.range.end; ++i)
                        {
                            for (int64_t j = i + 1; j < mappedReadsNum; ++j)
                            {
                                computeReadsPairCompatibility(result.statistics,
                                                              chunkAsqgFile,
                                                              chunkCsvEdgesFile,
                                                              i,
                                                              j,
                                                              bamRecordBuffer[i],
                                                              bamRecordBuffer[j],
                                                              readNames,
                                                              options);
                                ++result.completedPairs;

                                if (options.verbose && workerCount == 1)
                                {
                                    int64_t step = result.completedPairs <= result.range.pairCount / 100.0
                                        ? std::max<int64_t>(1, result.range.pairCount / 100000)
                                        : std::max<int64_t>(1, result.range.pairCount / 1000);
                                    printProgress(std::cout,
                                                  step,
                                                  result.completedPairs,
                                                  result.range.pairCount);
                                }
                            }
                        }

                        if (options.outputASQG)
                        {
                            chunkAsqgFile.close();
                            if (!chunkAsqgFile)
                                throw std::runtime_error("could not write ASQG chunk fragment");
                        }
                        if (options.outputCSV)
                        {
                            chunkCsvEdgesFile.close();
                            if (!chunkCsvEdgesFile)
                                throw std::runtime_error("could not write CSV chunk fragment");
                        }

                        result.statistics.numConsideredReadsPairs = result.completedPairs;
                    }
                    catch (std::exception const &error)
                    {
                        result.error = error.what();
                        stopWorkers.store(true);
                    }
                }
            });
        }
    }
    catch (...)
    {
        for (auto &worker : workers)
            worker.join();
        removeFragments(chunkResults);
        throw;
    }

    for (auto &worker : workers)
        worker.join();

    for (auto const &result : chunkResults)
    {
        if (!result.error.empty())
        {
            removeFragments(chunkResults);
            throw std::runtime_error(result.error);
        }
    }

    int64_t completedPairs = 0;
    for (auto const &result : chunkResults)
    {
        if (options.outputASQG && !appendFile(asqgFile, result.asqgFragment))
        {
            removeFragments(chunkResults);
            throw std::runtime_error("could not merge ASQG chunk fragment");
        }
        if (options.outputCSV && !appendFile(csvEdgesFile, result.csvEdgesFragment))
        {
            removeFragments(chunkResults);
            throw std::runtime_error("could not merge CSV chunk fragment");
        }

        completedPairs += result.completedPairs;
        mergeCompatibilityStatistics(globalStats, result.statistics);
    }

    removeFragments(chunkResults);
    asqgFile.close();
    csvEdgesFile.close();

    if (options.pairShardCount > 1)
    {
        std::string metadataFilename = outputBasename + ".pair-shard.tsv";
        std::ofstream metadataFile(metadataFilename,
                                   std::ofstream::out | std::ofstream::trunc);
        if (!metadataFile)
            throw std::runtime_error("could not open pair shard metadata file");

        metadataFile << "read_count\t" << mappedReadsNum << "\n"
                     << "total_pair_count\t" << maxReadsPairs << "\n"
                     << "shard_index\t" << options.pairShardIndex << "\n"
                     << "shard_count\t" << options.pairShardCount << "\n"
                     << "outer_begin\t" << shardRange.begin << "\n"
                     << "outer_end\t" << shardRange.end << "\n"
                     << "assigned_pair_count\t" << shardRange.pairCount << "\n"
                     << "completed_pair_count\t" << completedPairs << "\n"
                     << "requested_threads\t" << options.threads << "\n"
                     << "thread_limit\t" << MAX_WORKER_THREADS << "\n"
                     << "threads\t" << workerCount << "\n"
                     << "chunk_count\t" << chunkResults.size() << "\n";
    }

    if (options.verbose || options.pairShardCount > 1)
    {
        std::cout << "INFO: Pair shard " << options.pairShardIndex
                  << "/" << options.pairShardCount
                  << " assigned outer range [" << shardRange.begin
                  << ", " << shardRange.end << ")\n"
                  << "INFO: Assigned read pairs: " << shardRange.pairCount
                  << "/" << maxReadsPairs << "\n"
                  << "INFO: Completed read pairs: " << completedPairs << "\n\n";
    }
}

/******************************************************************************
    Initialize the graph vertices and store read names
******************************************************************************/
void initializeGraph(TGraph &graph,
                     std::vector<TVertexDescriptor > &vertices,
                     TProperties &readNames,
                     std::vector<std::vector<seqan::BamAlignmentRecord> > const &bamRecordBuffer,
                     AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    seqan::CharString oldQName = "";

    // Create a vertex for each read
    for (auto i=static_cast<unsigned>(0); i<bamRecordBuffer.size(); ++i)
    {
        vertices.push_back(addVertex(graph));
    }

    // Resize the property map containing the read names to fit all reads from the graph
    seqan::resizeVertexMap(readNames, graph);

    // Assign each read name to the corresponding vertex
    for (auto i=static_cast<unsigned>(0); i<bamRecordBuffer.size(); ++i)
    {
        seqan::assignProperty(readNames, vertices[i], bamRecordBuffer[i][0].qName);
    }
}

/******************************************************************************
    Write all reads sequences to ASQG output file
******************************************************************************/
void writeReadsToASQG(std::ostream &asqgFile,
                      std::vector<std::vector<seqan::BamAlignmentRecord> > const &bamRecordBuffer,
                      AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    int64_t i = 0;

    for (auto const &readBamRecordBuffer : bamRecordBuffer)
    {
        auto &readBamRecord = readBamRecordBuffer[0];

        seqan::CharString readSeq(readBamRecord.seq);

        if (seqan::hasFlagRC(readBamRecord))
        {
            seqan::reverseComplement(readSeq);
        }

        asqgFile << "VT\t" << i << "\t"
                 << readSeq << "\n";

        ++i;
    }
}

/******************************************************************************
    Write all reads sequences to CSV output file
******************************************************************************/
void writeReadsToCSV(std::ostream &csvNodesFile,
                     TProperties const &readNames,
                     std::vector<std::vector<seqan::BamAlignmentRecord> > const &bamRecordBuffer,
                     AlphaOptions const &options)
{
    if (options.debug)
    {
        std::cout << "DEBUG: Currently in file: " << __FILE__
                  << " Function: " << __FUNCTION__ << "()"
                  << "\n" << "\n" << std::flush;
    }

    csvNodesFile << "Id;Label;Specie" << "\n";

    int64_t i = 0;

    for (auto const &readBamRecordBuffer : bamRecordBuffer)
    {
        auto &readBamRecord = readBamRecordBuffer[0];

        auto specieID = seqan::prefix(readNames[i], 3);

        csvNodesFile << i << ";" << readBamRecord.qName << ";" << specieID << "\n";

        ++i;
    }
}

/******************************************************************************
    Compute the compatibility between 2 reads given all their bam records
******************************************************************************/
void computeReadsPairCompatibility(GlobalStatistics &globalStats,
                                   std::ofstream &asqgFile,
                                   std::ofstream &csvEdgesFile,
                                   int64_t i,
                                   int64_t j,
                                   std::vector<seqan::BamAlignmentRecord> const &readBamRecordBufferI,
                                   std::vector<seqan::BamAlignmentRecord> const &readBamRecordBufferJ,
                                   TProperties const &readNames,
                                   AlphaOptions const &options)
{
    bool atLeastOneCommonRef = false;
    bool areReadsCompatible = false;
    bool areReadsIncompatible = false;
    bool areReadsOverlapping = false;
    bool wasFoundWithMultiRef = false;

    bool trueCoEmittedReadsSpecie = (seqan::prefix(readNames[i], 3) == seqan::prefix(readNames[j], 3));

    // Initialize the overlap description
    OverlapDescription alignOvDescription;

    // Test each pair of bamRecords from read_i vs. read_j
    // and exit the double loop as soon as a compatible alignment
    // has been found

    // TO DO: Envisager de parcourir la double boucle en proposant les paires
    // d'alignements de bonne qualité en premier si possible
    // int32_t min_align_num = std::min(readBamRecordBufferI.size(), readBamRecordBufferJ.size());

    for (auto const &bamRecordI : readBamRecordBufferI)
    {
        for (auto const &bamRecordJ : readBamRecordBufferJ)
        {
            if (bamRecordI.rID == bamRecordJ.rID)
            {
                ++globalStats.numConsideredAlignmentsPairs;

                alignOvDescription.isRead2RC = !(seqan::hasFlagRC(bamRecordI)==seqan::hasFlagRC(bamRecordJ));

                atLeastOneCommonRef = true;

                computeAlignmentsPairCompatibility(areReadsCompatible,
                                                   areReadsIncompatible,
                                                   areReadsOverlapping,
                                                   globalStats,
                                                   bamRecordI,
                                                   bamRecordJ,
                                                   options);
            }

            // Exit the nested loop if there is enough information
            // to decide whether the reads are compatible or not
            if (areReadsCompatible || areReadsIncompatible)
            {
                goto fastExit;
            }
        }
    }

    // There we never went through the fastExit
    // these reads are neither compatible neither incompatible
    // TO DO: Understand what is happening. WTF ???
    // Si ni compatible ni incompatoble --> statut particulier
    // les mettre dans un paquet pour les traiter plus tard ?
    ++globalStats.neitherCompatNorIncompatReadPairsNum;
    //    std::cerr << "Reads are neither compatibles nor incompatibles..." << "\n";

    areReadsIncompatible = true;

// goto tag to exit the nested loop just above
fastExit:

    if (atLeastOneCommonRef)
        ++globalStats.readPairsWithCommonRefNum;

    // If at least one alignment pair was compatible
    // we add an edge between the 2 reads

    OverlapStatistics alignOvStats;

    if (areReadsCompatible)
    {
        if (areReadsIncompatible)
        {
            std::cerr << "ERROR: 2 reads cannot be both compatible and incompatible, go check your algo..."
                      << "\n" << "\n" << std::flush;
        }

        //
        auto const &bamRecordI = readBamRecordBufferI[0];
        auto const &bamRecordJ = readBamRecordBufferJ[0];

        // De novo overlap alignment computing
        auto readSeqI(bamRecordI.seq);
        if (seqan::hasFlagRC(bamRecordI)) seqan::reverseComplement(readSeqI);

        auto readSeqJ(bamRecordJ.seq);
        if (seqan::hasFlagRC(bamRecordJ) && alignOvDescription.isRead2RC)
        {
            //nothing
        }
        else if (seqan::hasFlagRC(bamRecordJ) || alignOvDescription.isRead2RC)
            seqan::reverseComplement(readSeqJ);

        // Perform the de-novo alignment of the 2 reads sequences, using the levenstein distance as a score scheme
        TAlign align;
        // int score = align2ReadSequences(align, readSeqI, readSeqJ);
        align2ReadSequences(align, readSeqI, readSeqJ);

        alignOvDescription.readIDI = i;
        alignOvDescription.readIDJ = j;

        computeAlignmentStats(alignOvStats, alignOvDescription, align, readSeqI, readSeqJ, false);

        auto seqIdentityPercent = static_cast<double>(alignOvStats.matchesNum) / alignOvStats.totalOverlapPositions;

        if (alignOvDescription.isEnclosedOverlap)
        {
            ++globalStats.enclosedOverlapsNum;
        }

        // Here are the major criteria to determine wether 2 reads are compatible.
        areReadsOverlapping = alignOvStats.totalOverlapPositions >= options.minOverlapLength;
        areReadsCompatible = (areReadsOverlapping && (alignOvStats.indelNum == 0)
                              && (seqIdentityPercent >= options.idRateThreshold));

        if (alignOvDescription.isRead2RC) reverseComplementReadJ(alignOvDescription);
    }

    if (areReadsCompatible)
    {
//        std::cerr << "Les 2 reads sont compatibles" << "\n";

        ++globalStats.compatibleReadPairsNum;

        // Add a new edge between the two reads
        if (options.outputASQG)
        {
            writeOverlapToASQG(asqgFile, alignOvDescription);
        }

        if (options.outputCSV)
        {
            writeOverlapToCSV(csvEdgesFile,
                              alignOvDescription,
                              trueCoEmittedReadsSpecie,
                              wasFoundWithMultiRef,
                              alignOvStats);
        }

        // Evaluate TP and FP
        if (trueCoEmittedReadsSpecie) ++globalStats.truePositive;
        else ++globalStats.falsePositive;
    }
    else // reads are considered incompatible by default
    {
        if (areReadsIncompatible)
        {
            ++globalStats.incompatibleReadPairsNum;
        }

//        std::cerr << "Les 2 reads ne sont pas compatibles" << "\n";
        // Evaluate FN et TN
        if(areReadsOverlapping)
        {
            if (trueCoEmittedReadsSpecie) ++globalStats.falseNegative;
            else ++globalStats.trueNegative;
        }
    }
}

/******************************************************************************
    Write an overlap to the ASQG output file
******************************************************************************/
void writeOverlapToASQG(std::ostream &asqgFile,
                        OverlapDescription const &aOvDcr)
{
    asqgFile << "ED\t" << aOvDcr.readIDI << " " << aOvDcr.readIDJ
             << " " << aOvDcr.ovBeginPosI << " " << aOvDcr.ovEndPosI
             << " " << aOvDcr.readLengthI
             << " " << aOvDcr.ovBeginPosJ << " " << aOvDcr.ovEndPosJ
             << " " << aOvDcr.readLengthJ << " " << aOvDcr.isRead2RC
             << " -1\n";
}

/******************************************************************************
    Write an overlap to the CSV Edges output file
******************************************************************************/
void writeOverlapToCSV(std::ostream &csvEdgesFile,
                       OverlapDescription const &aOvDcr,
                       bool const trueCoEmittedReadsSpecie,
                       bool const wasFoundWithMultiRef,
                       OverlapStatistics const &aOvStat)
{
    auto seqIdentityPercent = static_cast<double>(aOvStat.matchesNum) / aOvStat.totalOverlapPositions;

    csvEdgesFile << std::setprecision(4)
                 << aOvDcr.readIDI << ";" << aOvDcr.readIDJ << ";"
                 << "Undirected;"
                 << trueCoEmittedReadsSpecie << ";"
                 << (trueCoEmittedReadsSpecie ? 1 : 2) << ";"
                 << seqIdentityPercent << ";"
                 << wasFoundWithMultiRef << "\n";
}

/******************************************************************************
    Align 2 read sequences
******************************************************************************/
int align2ReadSequences(TAlign &align,
                        TSequence const &readSeqI,
                        TSequence const &readSeqJ)
{
    seqan::resize(seqan::rows(align), 2);
    seqan::assignSource(seqan::row(align, 0), readSeqI);
    seqan::assignSource(seqan::row(align, 1), readSeqJ);

//    auto scoringScheme = seqan::Score<int, seqan::Simple>(2, -4, -16, -3);
//    auto scoringScheme = seqan::Score<int, seqan::Simple>(1, -1, -2);
    auto scoringScheme = seqan::Score<int, seqan::Simple>(1, -1, -1);
    auto alignConfig = seqan::AlignConfig<true, true, true, true>();

//    int maxError = (int) (1 + std::min(length(readSeqI), length(readSeqJ)) * 0.05);

    int score = seqan::globalAlignment(align, scoringScheme, alignConfig);
//    int score = seqan::globalAlignment(align, scoringScheme, alignConfig, -maxError, maxError);

    return score;
}

/******************************************************************************
    Compute statistics for a pairwise alignment from Seqan
******************************************************************************/
void computeAlignmentStats(OverlapStatistics &alignOvStats,
                           OverlapDescription &alignOvDescription,
                           TAlign const &align,
                           TSequence const &readSeqI,
                           TSequence const &readSeqJ,
                           bool const toPrint)
{
    if (toPrint)
        std::cout << "\n" << align;

    auto &rowI = row(align, 0);
    auto &rowJ = row(align, 1);

    auto beginPosInAlignI = toViewPosition(rowI, 0);
    auto beginPosInAlignJ = toViewPosition(rowJ, 0);

    alignOvDescription.readLengthI = seqan::length(readSeqI);
    alignOvDescription.readLengthJ = seqan::length(readSeqJ);

    int lastPosI = seqan::length(readSeqI)-1;
    int lastPosJ = seqan::length(readSeqJ)-1;

    auto endPosInAlignI = toViewPosition(rowI, lastPosI);
    auto endPosInAlignJ = toViewPosition(rowJ, lastPosJ);

    int overlapBeginPos = std::max(beginPosInAlignI, beginPosInAlignJ);
    int overlapEndPos = std::min(endPosInAlignI, endPosInAlignJ);
    alignOvStats.totalOverlapPositions = overlapEndPos - overlapBeginPos + 1;

    alignOvDescription.ovBeginPosI = std::max(0, (int)toSourcePosition(rowI, beginPosInAlignJ));
    alignOvDescription.ovBeginPosJ = std::max(0, (int)toSourcePosition(rowJ, beginPosInAlignI));

    alignOvDescription.ovEndPosI = std::min(lastPosI, (int)toSourcePosition(rowI, endPosInAlignJ));
    alignOvDescription.ovEndPosJ = std::min(lastPosJ, (int)toSourcePosition(rowJ, endPosInAlignI));

    if ((alignOvDescription.ovBeginPosI == 0 && alignOvDescription.ovEndPosI == lastPosI)
        || (alignOvDescription.ovBeginPosJ == 0 && alignOvDescription.ovEndPosJ == lastPosJ))
    {
        alignOvDescription.isEnclosedOverlap = true;
    }

    if (!alignOvDescription.isEnclosedOverlap)
    {
        for (int i = overlapBeginPos; i <= overlapEndPos; ++i)
        {
            if (seqan::isGap(rowI, i) || seqan::isGap(rowJ, i))
            {
                ++alignOvStats.indelNum;
            }
            else
            {
                if (rowI[i] == rowJ[i])
                {
                    ++alignOvStats.matchesNum;
                }
                else
                {
                    ++alignOvStats.mismatchesNum;
                }
            }

            if (toPrint)
            {
                std::cout << i << " " << rowI[i] << " " << rowJ[i] << " "
                          << alignOvStats.matchesNum << " "
                          << alignOvStats.mismatchesNum << " "
                          << alignOvStats.indelNum << "\n";
            }
        }
    }

    if (toPrint)
    {
        std::cout << "\n" << alignOvStats.totalOverlapPositions << "\n"
                  << alignOvStats.matchesNum << "\n"
                  << alignOvStats.indelNum << "\n";
    }
}

/******************************************************************************
    Print the progress of the nested loop
******************************************************************************/
void printProgress(std::ostream &stream,
                   int64_t step,
                   int64_t numConsideredReadsPairs,
                   int64_t maxReadsPairs)
{
    if (numConsideredReadsPairs % step == 0)
    {
        double progress = (double) numConsideredReadsPairs * 100.0 / maxReadsPairs;
        stream << std::setprecision(5)
               << "\rINFO: Compatibility graph building progress at "
               << progress << "%         " << std::flush;
    }
}
