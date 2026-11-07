#!/usr/bin/env python3

import gzip

def read_fastq_file_handle(fastq_file_handle):
    """
    Parse a fastq file and return a generator
    """
    # Variables initialization
    count = 0
    header = ''
    seq = ''
    qual = ''
    # Reading input file
    for line in (l.strip() for l in fastq_file_handle if l.strip()):
        count += 1
        if count % 4 == 1:
            if header:
                yield header, seq, qual
            header = line[1:].split()[0]
        elif count % 4 == 2:
            seq = line
        elif count % 4 == 0:
            qual = line
    # yield last fastq sequence (nothing for an empty file)
    if header:
        yield header, seq, qual
    # Close input file
    fastq_file_handle.close()


# Lowest quality character of the Phred+64 (Illumina 1.3-1.7) and Solexa+64
# encodings: Solexa scores go down to -5, i.e. ';' (59).
PHRED64_MIN_ASCII = 59
# Highest character of a typical Illumina Phred+33 file ('J', Q41). Above it
# both encodings are possible (Phred+33 Q42+, Phred+64 Q11+).
PHRED33_ILLUMINA_MAX_ASCII = 74


def is_phred33(fastq, number_of_reads_to_test=-1):
    """
    Determine if the fastq uses an offset of +33.
    See for more details: https://en.wikipedia.org/wiki/FASTQ_format

    A quality character below ';' (59) can only be Phred+33, so it settles
    the question. Without one, the file is still accepted as Phred+33 when
    every character is within the usual Illumina Phred+33 range ('!'..'J'),
    e.g. uniformly high qualities ('IIII'). Only a file with no character
    below ';' and some above 'J', the profile of Phred+64 data, is rejected,
    as is any character outside the printable FASTQ range ('!'..'~').
    An empty file is not rejected.
    By default, test until end of file.
    """
    min_ascii_code = None
    max_ascii_code = None
    opener = gzip.open if fastq.endswith('.gz') else open
    with opener(fastq, 'rt') as handler:
        for i, (header, seq, qual) in enumerate(read_fastq_file_handle(handler)):
            if number_of_reads_to_test > 0 and i >= number_of_reads_to_test:
                break
            if not qual:
                continue
            ascii_codes = [ord(qual_char) for qual_char in qual]
            current_min, current_max = min(ascii_codes), max(ascii_codes)
            if current_min < 33 or current_max > 126:
                return False
            if current_min < PHRED64_MIN_ASCII:
                return True
            min_ascii_code = current_min if min_ascii_code is None else min(current_min, min_ascii_code)
            max_ascii_code = current_max if max_ascii_code is None else max(current_max, max_ascii_code)

    if max_ascii_code is None:
        return True
    return max_ascii_code <= PHRED33_ILLUMINA_MAX_ASCII
