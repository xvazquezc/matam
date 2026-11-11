#!/usr/bin/env python3

import gzip
import itertools
import queue
import threading
import zlib

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


def _phred33_verdict(min_ascii_code, max_ascii_code):
    """True if qualities spanning these ASCII codes can be Phred+33 (see is_phred33)."""
    if min_ascii_code < 33 or max_ascii_code > 126:
        return False
    return min_ascii_code < PHRED64_MIN_ASCII or max_ascii_code <= PHRED33_ILLUMINA_MAX_ASCII


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


FASTQ_SUFFIXES = ('.fq', '.fastq', '.fq.gz', '.fastq.gz')
_MATE_SEPARATORS = '/._:'


class FastqFormatError(ValueError):
    """An input FASTQ file is malformed, or the R1/R2 files don't hold the same pairs in order."""


def is_same_pair(id1, id2):
    """True if R1 and R2 read IDs (first header tokens) name the same pair.

    Mates share the ID (CASAVA 1.8+, SRA), or the IDs differ only by a trailing
    mate number 1/2 after the same separator: ``r/1`` ``r/2`` (legacy),
    ``SRR1.7.1`` ``SRR1.7.2`` (fastq-dump --readids), ``r_1`` ``r_2``.
    """
    return id1 == id2 or (len(id1) == len(id2) > 2 and id1[:-1] == id2[:-1]
                          and id1[-2] in _MATE_SEPARATORS and id1[-1] == '1' and id2[-1] == '2')


def _text(value):
    return value.decode('utf-8', 'replace') if isinstance(value, bytes) else value


def _check_record(lines, path, number):
    """Check one record (4 byte lines, with their line ends); return its read ID."""
    header, sequence, plus, quality = lines
    if not header.strip():
        raise FastqFormatError('{}: blank line at record {}'.format(path, number))
    if header[:1] != b'@' or plus[:1] != b'+':
        raise FastqFormatError(
            '{}: record {} is not FASTQ (expected an "@" header line and a "+" line): '
            '{!r}'.format(path, number, _text(header.rstrip(b'\r\n'))))
    sequence, quality = sequence.rstrip(b'\r\n'), quality.rstrip(b'\r\n')
    if len(sequence) != len(quality):
        raise FastqFormatError('{}: record {} has {} bases but {} quality values'.format(
            path, number, len(sequence), len(quality)))
    fields = header[1:].split(None, 1)
    return fields[0] if fields else b''


_BATCH = 10000   # records checked together
# Pairs whose qualities are checked: the encoding is the same for a whole
# file, so the first reads tell it (the old check used 400 reads per file).
_QUALITY_SAMPLE = 1000
# Quality bytes: printable FASTQ range, Phred+33-only range, and above 'J'.
_LINE_ENDS = b'\r\n'
_VALID_QUALITY = bytes(range(33, 127)) + _LINE_ENDS
_PHRED33_ONLY = bytes(range(33, PHRED64_MIN_ASCII))
_ABOVE_ILLUMINA_PHRED33 = bytes(range(PHRED33_ILLUMINA_MAX_ASCII + 1, 127))


def _open_binary(path):
    """A binary stream of a plain or gzipped file.

    zlib releases the GIL while it decompresses, so reading each file in its
    own thread (:class:`_Prefetch`) decompresses the two files of a pair in
    parallel with the checks.
    """
    if path.lower().endswith('.gz'):
        return gzip.open(path, 'rb')
    return open(path, 'rb')


def _batches(handle, path):
    """Yield the lines of a FASTQ stream in batches of whole records (4 lines each).

    Blank lines are allowed at the end of the file only; a file ending inside a
    record raises FastqFormatError.
    """
    records = 0
    while True:
        try:
            lines = list(itertools.islice(handle, 4 * _BATCH))
        except (OSError, EOFError, zlib.error) as error:   # e.g. a truncated .gz file
            raise FastqFormatError('{}: cannot be read: {}'.format(path, error))
        if len(lines) < 4 * _BATCH:          # the end of the file
            while lines and not lines[-1].strip():
                lines.pop()
            if lines and not lines[-1].endswith(b'\n'):
                lines[-1] += b'\n'
            if len(lines) % 4:
                blank = next((i for i, line in enumerate(lines) if not line.strip()), None)
                if blank is not None:
                    raise FastqFormatError('{}: blank line at record {}'.format(
                        path, records + blank // 4 + 1))
                raise FastqFormatError(
                    '{}: the last record is truncated (the file ends inside it)'.format(path))
            if lines:
                yield lines
            return
        records += _BATCH
        yield lines


class _Prefetch:
    """Iterate over an iterator in a background thread, a few items ahead.

    Reading a pipe or a file releases the GIL, so the decompressor of each file
    keeps running while the other file is checked.
    """
    _END = object()

    def __init__(self, iterator, ahead=4):
        self._queue = queue.Queue(maxsize=ahead)
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, args=(iterator,), daemon=True)
        self._thread.start()

    def _put(self, item):
        while not self._stop.is_set():
            try:
                self._queue.put(item, timeout=0.1)
                return True
            except queue.Full:
                pass
        return False

    def _run(self, iterator):
        try:
            for item in iterator:
                if not self._put((True, item)):
                    return
        except BaseException as error:      # handed to the consumer
            self._put((False, error))
            return
        self._put((True, self._END))

    def __iter__(self):
        while True:
            ok, item = self._queue.get()
            if not ok:
                raise item
            if item is self._END:
                return
            yield item

    def close(self):
        self._stop.set()
        self._thread.join()


def _batch_is_well_formed(lines):
    """True if every record of a batch has '@', '+' and as many quality values as bases."""
    return (all(map(bytes.startswith, lines[0::4], itertools.repeat(b'@')))
            and all(map(bytes.startswith, lines[2::4], itertools.repeat(b'+')))
            and list(map(len, lines[1::4])) == list(map(len, lines[3::4])))


def _read_ids(headers):
    return [fields[0] if fields else b'' for fields in
            map(bytes.split, [header[1:] for header in headers], itertools.repeat(None),
                itertools.repeat(1))]


def _mates_by_suffix(ids1, ids2):
    """True if every R1/R2 ID pair differs only by a trailing 1/2 after a separator."""
    return (all(len(i) > 2 for i in ids1)
            and [i[:-1] for i in ids1] == [i[:-1] for i in ids2]
            and all(map(bytes.endswith, ids1, itertools.repeat(b'1')))
            and all(map(bytes.endswith, ids2, itertools.repeat(b'2')))
            and {i[-2:-1] for i in ids1} <= {c.encode() for c in _MATE_SEPARATORS})


def validate_paired_fastq(forward, reverse):
    """Check the R1/R2 FASTQ input of a run; return its number of read pairs.

    Both files are read once, in step, and every record is checked: FASTQ
    structure (``@`` header, ``+`` line, as many quality values as bases, no
    truncated last record), the same number of records in both files, mates in
    the same order (:func:`is_same_pair` on every pair: MATAM pairs mates by
    position). The quality encoding, the same for a whole file, is checked on
    the first ``_QUALITY_SAMPLE`` pairs: Phred+33 by the rule of
    :func:`is_phred33` (SGA needs Phred+33). Plain or gzipped (``.gz``) files,
    each read (and decompressed) in its own thread; records are checked in
    batches with C-level operations. Raises FastqFormatError naming the file,
    the record and the problem.
    """
    for path in (forward, reverse):
        if not path.lower().endswith(FASTQ_SUFFIXES):
            raise FastqFormatError('{}: not a FASTQ file name (expected {})'.format(
                path, ', '.join(FASTQ_SUFFIXES)))
    with _open_binary(forward) as forward_handle, _open_binary(reverse) as reverse_handle:
        forward_batches = _Prefetch(_batches(forward_handle, forward))
        reverse_batches = _Prefetch(_batches(reverse_handle, reverse))
        try:
            pairs, phred33_only, above_phred33 = _check_batches(
                forward, reverse, forward_batches, reverse_batches)
        finally:
            forward_batches.close()
            reverse_batches.close()
    if pairs == 0:
        raise FastqFormatError('{} and {} contain no reads'.format(forward, reverse))
    if above_phred33 and not phred33_only:
        raise FastqFormatError(
            '{} / {}: no quality value is below {!r} and some are above {!r}: this looks like '
            'Phred+64, not the Phred+33 that SGA needs'.format(
                forward, reverse, chr(PHRED64_MIN_ASCII), chr(PHRED33_ILLUMINA_MAX_ASCII)))
    return pairs


def _check_batches(forward, reverse, forward_batches, reverse_batches):
    """The checks of validate_paired_fastq; returns (pairs, phred33_only, above_phred33)."""
    pairs = 0
    phred33_only = above_phred33 = False
    for batch1, batch2 in itertools.zip_longest(forward_batches, reverse_batches):
        if batch1 is None or batch2 is None or len(batch1) != len(batch2):
            records1 = pairs + (len(batch1) // 4 if batch1 else 0)
            records2 = pairs + (len(batch2) // 4 if batch2 else 0)
            longer = forward if records1 > records2 else reverse
            raise FastqFormatError(
                'The paired FASTQ files have different numbers of reads: {} has more than '
                'the {} reads of the other file'.format(longer, min(records1, records2)))
        for batch, path in ((batch1, forward), (batch2, reverse)):
            if not _batch_is_well_formed(batch):
                # Find and report the faulty record.
                for index in range(0, len(batch), 4):
                    _check_record(batch[index:index + 4], path, pairs + index // 4 + 1)
        ids1, ids2 = _read_ids(batch1[0::4]), _read_ids(batch2[0::4])
        if ids1 != ids2 and not _mates_by_suffix(ids1, ids2):
            for index, (id1, id2) in enumerate(zip(ids1, ids2)):
                if not is_same_pair(_text(id1), _text(id2)):
                    raise FastqFormatError(
                        'The paired FASTQ files are not in step: read {} is {!r} in {} but '
                        '{!r} in {}. MATAM pairs mates by position, so both files must list '
                        'the same pairs in the same order'.format(
                            pairs + index + 1, _text(id1), forward, _text(id2), reverse))
        for batch, path in ((batch1, forward), (batch2, reverse)):
            if pairs >= _QUALITY_SAMPLE:
                break
            qualities = b''.join(batch[3::4][:_QUALITY_SAMPLE - pairs])
            if qualities.translate(None, _VALID_QUALITY):
                for index in range(3, len(batch), 4):
                    bad = batch[index].translate(None, _VALID_QUALITY)
                    if bad:
                        raise FastqFormatError(
                            '{}: record {} has a quality character outside the FASTQ '
                            'range: {!r}'.format(path, pairs + index // 4 + 1, _text(bad[:1])))
            if not phred33_only:
                phred33_only = len(qualities.translate(None, _PHRED33_ONLY)) != len(qualities)
                if not phred33_only and not above_phred33:
                    above_phred33 = len(qualities.translate(None, _ABOVE_ILLUMINA_PHRED33)) != len(qualities)
        pairs += len(batch1) // 4
    return pairs, phred33_only, above_phred33
