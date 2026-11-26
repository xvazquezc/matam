#!/usr/bin/env python3
"""SAM helpers shared by the MATAM scripts."""


def parse_cigar(cigar):
    """
    Parse a CIGAR string and return a list of (operation, count) tuples.
    A missing count means 1, and adjacent repeats of an operation are merged
    ('5M5M' -> [('M', 10)]).
    """
    cigar_tab = list()
    count_str = ''
    for c in cigar:
        if c.isdigit():
            count_str += c
        else:
            operation = c
            count = 1
            if count_str:
                count = int(count_str)
            if cigar_tab and operation == cigar_tab[-1][0]:
                cigar_tab[-1] = (operation, cigar_tab[-1][1] + count)
            else:
                cigar_tab.append((operation, count))
            count_str = ''
    return cigar_tab
