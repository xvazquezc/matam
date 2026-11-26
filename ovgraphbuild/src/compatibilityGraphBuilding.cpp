#include "compatibilityGraphBuilding.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>

namespace
{
constexpr int MAX_WORKER_THREADS = 24;
constexpr int CHUNKS_PER_WORKER = 8;

/******************************************************************************
    Which reads share a reference.

    Two reads are compared alignment by alignment, and only alignments on the
    same reference can make them compatible or incompatible. A pair of reads
    with no common reference is always "neither compatible nor incompatible"
    (then counted as incompatible) and never gives an edge. Indexing the reads
    by reference lets the pair loop visit only the pairs with a common
    reference and account for the others in bulk.
******************************************************************************/
struct ReferenceIndex
{
    // Read indices aligned to each reference, ascending. Bucket rID + 1, so
    // that a record without a reference (rID -1) keeps its own bucket.
    std::vector<std::vector<int32_t> > readsByRefBucket;
    // Distinct reference buckets of each read.
    std::vector<std::vector<int32_t> > refBucketsByRead;
};

ReferenceIndex buildReferenceIndex(std::vector<std::vector<seqan::BamAlignmentRecord> > const &bamRecordBuffer)
{
    if (bamRecordBuffer.size() > static_cast<size_t>(std::numeric_limits<int32_t>::max()))
        throw std::runtime_error("too many mapped reads to index by reference");

    ReferenceIndex index;
    index.refBucketsByRead.resize(bamRecordBuffer.size());

    int32_t maxBucket = 0;
    for (size_t read = 0; read < bamRecordBuffer.size(); ++read)
    {
        auto &buckets = index.refBucketsByRead[read];
        for (auto const &bamRecord : bamRecordBuffer[read])
        {
            if (bamRecord.rID < -1)
                throw std::runtime_error("invalid reference ID in SAM record");
            buckets.push_back(bamRecord.rID + 1);
        }
        std::sort(buckets.begin(), buckets.end());
        buckets.erase(std::unique(buckets.begin(), buckets.end()), buckets.end());
        if (!buckets.empty())
            maxBucket = std::max(maxBucket, buckets.back());
    }

    index.readsByRefBucket.resize(static_cast<size_t>(maxBucket) + 1);
    for (size_t read = 0; read < bamRecordBuffer.size(); ++read)
        for (auto bucket : index.refBucketsByRead[read])
            index.readsByRefBucket[bucket].push_back(static_cast<int32_t>(read));

    return index;
}

/******************************************************************************
    Reads j > i that share at least one reference with read i, ascending.
    lastSeen[j] == i marks j as already collected for this i.
******************************************************************************/
void collectCandidatePartners(std::vector<int32_t> &candidates,
                              std::vector<int32_t> &lastSeen,
                              ReferenceIndex const &index,
                              int32_t i)
{
    candidates.clear();
    for (auto bucket : index.refBucketsByRead[i])
    {
        auto const &reads = index.readsByRefBucket[bucket];
        for (auto it = std::upper_bound(reads.begin(), reads.end(), i); it != reads.end(); ++it)
        {
            if (lastSeen[*it] != i)
            {
                lastSeen[*it] = i;
                candidates.push_back(*it);
            }
        }
    }
    std::sort(candidates.begin(), candidates.end());
}

struct PairRange
{
    int64_t begin;
    int64_t end;
    int64_t pairCount;
    int64_t work;
};

struct PairChunkResult
{
    PairRange range;
    GlobalStatistics statistics;
    int64_t completedPairs = 0;
    int64_t comparedPairs = 0;
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

/******************************************************************************
    Expected work of the outer indices, to balance shards and thread chunks.

    workBefore[k] is the work of the outer indices 0..k-1 (k = 0..n-1); a
    range [begin, end) costs workBefore[end] - workBefore[begin].

    - All pairs: read i is compared with its n-1-i partners, so the work is
      the pair count.
    - Common-reference pairs (default): read i costs one alignment-pair
      comparison per (later read, shared reference), plus 1 for the read
      itself (collecting the partners, accounting for the skipped pairs).
      Balancing by the pair count instead put most of the comparisons in a
      few ranges: reads differ widely in how many reads share a reference.
******************************************************************************/
std::vector<int64_t> outerWorkPrefix(int64_t readCount,
                                     bool allPairs,
                                     ReferenceIndex const &index)
{
    int64_t outerCount = std::max<int64_t>(0, readCount - 1);
    std::vector<int64_t> workBefore(static_cast<size_t>(outerCount) + 1, 0);
    for (int64_t i = 0; i < outerCount; ++i)
    {
        int64_t work;
        if (allPairs)
        {
            work = readCount - 1 - i;
        }
        else
        {
            work = 1;
            for (auto bucket : index.refBucketsByRead[i])
            {
                auto const &reads = index.readsByRefBucket[bucket];
                work += reads.end() - std::upper_bound(reads.begin(), reads.end(),
                                                       static_cast<int32_t>(i));
            }
        }
        workBefore[i + 1] = workBefore[i] + work;
    }
    return workBefore;
}

int64_t nearestOuterBoundary(std::vector<int64_t> const &workBefore,
                             int64_t minOuterIndex,
                             int64_t maxOuterIndex,
                             int64_t targetWork)
{
    int64_t low = minOuterIndex;
    int64_t high = maxOuterIndex;

    while (low < high)
    {
        int64_t mid = low + (high - low) / 2;
        if (workBefore[mid] < targetWork)
            low = mid + 1;
        else
            high = mid;
    }

    if (low > minOuterIndex)
    {
        int64_t previous = low - 1;
        int64_t previousDistance = targetWork - workBefore[previous];
        int64_t currentDistance = workBefore[low] - targetWork;
        if (previousDistance <= currentDistance)
            return previous;
    }

    return low;
}

PairRange partitionPairRange(int64_t readCount,
                             std::vector<int64_t> const &workBefore,
                             int64_t outerBegin,
                             int64_t outerEnd,
                             int partIndex,
                             int partCount)
{
    int64_t workBegin = workBefore[outerBegin];
    int64_t workEnd = workBefore[outerEnd];
    int64_t beginTarget = balancedTarget(workBegin, workEnd, partIndex, partCount);
    int64_t endTarget = balancedTarget(workBegin, workEnd, partIndex + 1, partCount);

    int64_t begin = nearestOuterBoundary(workBefore, outerBegin, outerEnd, beginTarget);
    int64_t end = nearestOuterBoundary(workBefore, begin, outerEnd, endTarget);

    return PairRange{begin, end, pairsBefore(readCount, end) - pairsBefore(readCount, begin),
                     workBefore[end] - workBefore[begin]};
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
    destination.realignedReadPairsNum += source.realignedReadPairsNum;
    destination.perfectOverlapShortcutsNum += source.perfectOverlapShortcutsNum;
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
    ReferenceIndex referenceIndex;
    if (!options.allPairs)
        referenceIndex = buildReferenceIndex(bamRecordBuffer);
    std::vector<int64_t> workBefore = outerWorkPrefix(mappedReadsNum,
                                                      options.allPairs,
                                                      referenceIndex);
    PairRange shardRange = partitionPairRange(mappedReadsNum,
                                              workBefore,
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
                                             workBefore,
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
                std::vector<int32_t> candidates;
                std::vector<int32_t> lastSeen;
                if (!options.allPairs)
                    lastSeen.assign(mappedReadsNum, -1);

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

                        int64_t progressStep = std::max<int64_t>(1, result.range.pairCount / 1000);
                        int64_t nextProgress = 0;

                        for (int64_t i = result.range.begin; i < result.range.end; ++i)
                        {
                            int64_t partnersNum = mappedReadsNum - 1 - i;

                            if (options.allPairs)
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
                                }
                                result.comparedPairs += partnersNum;
                            }
                            else
                            {
                                collectCandidatePartners(candidates,
                                                         lastSeen,
                                                         referenceIndex,
                                                         static_cast<int32_t>(i));
                                for (auto j : candidates)
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
                                }

                                // The skipped pairs have no common reference:
                                // each one is counted exactly as the
                                // comparison would count it.
                                int64_t skippedNum = partnersNum - static_cast<int64_t>(candidates.size());
                                result.statistics.neitherCompatNorIncompatReadPairsNum += skippedNum;
                                result.statistics.incompatibleReadPairsNum += skippedNum;
                                result.comparedPairs += candidates.size();
                            }
                            result.completedPairs += partnersNum;

                            if (options.verbose && workerCount == 1
                                && result.completedPairs >= nextProgress)
                            {
                                printProgress(std::cout,
                                              1,
                                              result.completedPairs,
                                              result.range.pairCount);
                                nextProgress = result.completedPairs + progressStep;
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
    int64_t comparedPairs = 0;
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
        comparedPairs += result.comparedPairs;
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
                     << "assigned_work\t" << shardRange.work << "\n"
                     << "total_work\t" << workBefore[outerEnd] << "\n"
                     << "completed_pair_count\t" << completedPairs << "\n"
                     << "compared_pair_count\t" << comparedPairs << "\n"
                     << "all_pairs\t" << (options.allPairs ? 1 : 0) << "\n"
                     << "full_realignment\t" << (options.fullRealignment ? 1 : 0) << "\n"
                     << "any_shared_ref\t" << (options.anySharedRef ? 1 : 0) << "\n"
                     << "realigned_pair_count\t" << globalStats.realignedReadPairsNum << "\n"
                     << "perfect_overlap_shortcut_count\t" << globalStats.perfectOverlapShortcutsNum << "\n"
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
                  << "INFO: Completed read pairs: " << completedPairs << "\n"
                  << "INFO: Compared read pairs: " << comparedPairs
                  << (options.allPairs ? " (all pairs)" : " (with a common reference)")
                  << "\n\n";
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

    // A node is one read: (Label, Mate) is unique. Mate is 1/2 for paired
    // reads and 0 for unpaired ones. ComponentSearch copies these lines through.
    csvNodesFile << "Id;Label;Mate;Specie" << "\n";

    int64_t i = 0;

    for (auto const &readBamRecordBuffer : bamRecordBuffer)
    {
        auto &readBamRecord = readBamRecordBuffer[0];

        auto specieID = seqan::prefix(readNames[i], 3);

        csvNodesFile << i << ";" << readBamRecord.qName << ";" << mateOf(readBamRecord)
                     << ";" << specieID << "\n";

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
            // to decide whether the reads are compatible or not. By default
            // the first shared reference decides; with --any-shared-ref only
            // a compatible alignment pair ends the search.
            if (areReadsCompatible)
            {
                areReadsIncompatible = false;
                goto fastExit;
            }
            if (areReadsIncompatible && !options.anySharedRef)
            {
                goto fastExit;
            }
        }
    }

    // --any-shared-ref: every shared reference was incompatible
    if (areReadsIncompatible)
    {
        goto fastExit;
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

        alignOvDescription.readIDI = i;
        alignOvDescription.readIDJ = j;

        ++globalStats.realignedReadPairsNum;
        OverlapDescription const descriptionBefore = alignOvDescription;
        auto realignWithSeqAn = [&readSeqI, &readSeqJ](OverlapStatistics &stats,
                                                       OverlapDescription &description)
        {
            // Perform the de-novo alignment of the 2 reads sequences, using the levenstein distance as a score scheme
            TAlign align;
            align2ReadSequences(align, readSeqI, readSeqJ);
            computeAlignmentStats(stats, description, align, readSeqI, readSeqJ, false);
        };

        // Reads overlapping without any difference get the same statistics
        // without the full alignment (see perfectOverlap.h).
        if (!options.fullRealignment
            && computePerfectOverlap(alignOvStats, alignOvDescription, readSeqI, readSeqJ))
        {
            ++globalStats.perfectOverlapShortcutsNum;
            if (options.checkPerfectOverlaps)
            {
                OverlapStatistics seqanStats;
                OverlapDescription seqanDescription = descriptionBefore;
                realignWithSeqAn(seqanStats, seqanDescription);
                if (!sameOverlap(alignOvStats, alignOvDescription, seqanStats, seqanDescription))
                {
                    throw std::runtime_error("perfect-overlap shortcut differs from the SeqAn alignment for reads "
                                             + std::to_string(i) + " and " + std::to_string(j));
                }
            }
        }
        else
        {
            realignWithSeqAn(alignOvStats, alignOvDescription);
        }

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
