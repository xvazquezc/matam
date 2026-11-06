#ifndef SEQANOBJECTSFORMATING_H
#define SEQANOBJECTSFORMATING_H

#include <stdlib.h>
#include <stdint.h>

#include <seqan/bam_io.h>

/***************************************************************************************
    Format a bamRecord and return it as a string
***************************************************************************************/
std::string formatBamRecordString(seqan::BamAlignmentRecord const &bamRecord);

/***************************************************************************************
    Mate of a paired record: 1 (flag 0x40), 2 (flag 0x80), or 0 when unpaired.
    MATAM sets these flags after mapping, so that mates sharing a QNAME stay
    distinct reads.
***************************************************************************************/
inline int mateOf(seqan::BamAlignmentRecord const &bamRecord)
{
    bool const first = seqan::hasFlagFirst(bamRecord);
    bool const last = seqan::hasFlagLast(bamRecord);
    if (first == last)
        return 0;
    return first ? 1 : 2;
}

#endif // SEQANOBJECTSFORMATING_H
