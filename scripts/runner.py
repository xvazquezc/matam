#!/usr/bin/env python3

import sys
import os
import logging
import subprocess

# get the root handler to update his behavior (No prefix, no endline)
logger = logging.getLogger()


def logged_call(command, verbose=False, pipefail=False):
    """
    A logged version of subprocess.call. Do not wait the end of the
    process to start logging.

    With pipefail=True the command runs under bash with ``set -o pipefail``,
    so a failure anywhere in a shell pipeline is reported, not only a
    failure of its last command.
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
        null = open(os.devnull, 'w')
        stdout = null
        stderr = null

    #Update temporally the handler (No prefix, no endline)
    old_handlers = {handler: handler.formatter for handler in logger.handlers}
    for handler in logger.handlers:
        handler.formatter = logging.Formatter('%(message)s')
        handler.terminator = ''

    with subprocess.Popen(command, stdout=stdout, stderr=stderr, shell=True, bufsize=0,
                          executable=executable) as process:
        if verbose:
            while process.poll() is None:
                message = os.read(process.stdout.fileno(), 1024).decode()
                logger.info(message)


    # rehabilitate previous handler
    for handler in logger.handlers:
        handler.formatter = old_handlers[handler]
        handler.terminator = '\n'

    return process.wait()


def logged_check_call(command, verbose=False, pipefail=False):
    """
    A logged version of subprocess.check_call
    """
    returncode = logged_call(command, verbose=verbose, pipefail=pipefail)
    if returncode != 0:
        logger.info('CMD: {0}'.format(command))
        logger.fatal('The last command returns a non-zero return code: %s' % returncode)
        sys.exit('Non-zero return code')
