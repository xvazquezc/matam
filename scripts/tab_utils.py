#!/usr/bin/env python3
"""Reading of sorted, whitespace-separated tables shared by the MATAM scripts."""


def group_rows(rows, factor_index=0):
    """
    Yield the runs of consecutive rows (lists of fields) that share the value
    of column factor_index, each run as a list. No rows yield nothing.
    """
    previous_factor = None
    factor_tab_list = list()
    for tab in rows:
        current_factor = tab[factor_index]
        if current_factor != previous_factor and factor_tab_list:
            yield factor_tab_list
            factor_tab_list = list()
        factor_tab_list.append(tab)
        previous_factor = current_factor
    if factor_tab_list:
        yield factor_tab_list


def read_tab_file_handle_sorted(tab_file_handle, factor_index=0):
    """
    Parse a whitespace-separated table sorted by one column, and yield the
    rows of each run of equal values in that column as a list of split rows
    (see group_rows).

    Blank lines are skipped and an empty input yields nothing. The input can
    be a file handle or any iterable of lines; it is closed at the end when it
    has a close() method.
    """
    rows = (tab for tab in (l.split() for l in tab_file_handle) if tab)
    yield from group_rows(rows, factor_index)
    close = getattr(tab_file_handle, 'close', None)
    if close is not None:
        close()
