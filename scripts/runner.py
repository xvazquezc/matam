#!/usr/bin/env python3

import codecs
import sys
import os
import logging
import subprocess

# get the root handler to update his behavior (No prefix, no endline)
logger = logging.getLogger()


def logged_call(command, verbose=False, pipefail=True):
    """
    A logged version of subprocess.call. Do not wait the end of the
    process to start logging.

    By default the command runs under bash with ``set -o pipefail``, so a
    failure anywhere in a shell pipeline is reported, not only a failure of
    its last command. A pipeline whose inner commands may legitimately exit
    non-zero (e.g. ``grep`` with no match) must pass pipefail=False.
    """
    logger.debug('CMD: {0}'.format(command))
    executable = None
    if pipefail:
        command = 'set -o pipefail; ' + command
        executable = '/bin/bash'

    # if true, redirect stderr to stdout, otherwise, to NULL
    if verbose:
        stdout = subprocess.PIPE
        stderr = subprocess.STDOUT
    else:
        stdout = subprocess.DEVNULL
        stderr = subprocess.DEVNULL

    #Update temporally the handler (No prefix, no endline)
    old_handlers = {handler: handler.formatter for handler in logger.handlers}
    for handler in logger.handlers:
        handler.formatter = logging.Formatter('%(message)s')
        handler.terminator = ''

    with subprocess.Popen(command, stdout=stdout, stderr=stderr, shell=True, bufsize=0,
                          executable=executable) as process:
        if verbose:
            # Read until EOF, not until the process exits, so the end of the
            # output is not lost; decode incrementally so a multi-byte
            # character split across two reads stays intact.
            decoder = codecs.getincrementaldecoder('utf-8')(errors='replace')
            for chunk in iter(lambda: os.read(process.stdout.fileno(), 1024), b''):
                logger.info(decoder.decode(chunk))
            tail = decoder.decode(b'', final=True)
            if tail:
                logger.info(tail)


    # rehabilitate previous handler
    for handler in logger.handlers:
        handler.formatter = old_handlers[handler]
        handler.terminator = '\n'

    return process.wait()


def logged_check_call(command, verbose=False, pipefail=True):
    """
    A logged version of subprocess.check_call
    """
    returncode = logged_call(command, verbose=verbose, pipefail=pipefail)
    if returncode != 0:
        logger.info('CMD: {0}'.format(command))
        logger.fatal('The last command returns a non-zero return code: %s' % returncode)
        sys.exit('Non-zero return code')
