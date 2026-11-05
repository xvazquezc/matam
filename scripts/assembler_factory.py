#!/usr/bin/env python3

import binary_utils
import os
import sys
import logging
import runner
import shutil

logger = logging.getLogger(__name__)


class Assembler:

    @classmethod
    def name(cls):
        return cls.__name__


    def __init__(self):
        self.assembler_wrapper = self._assembler_wrapper()
        self.assembler_bin = self._assembler_bin()
        self.cmd_line = None

    def _assembler_wrapper(self):
        """
        return the assembler wrapper path
        """
        raise NotImplementedError( "Should have implemented this method: '%s'" % self._assembler_wrapper.__name__ )


    def _assembler_bin(self):
        """
        return the assembler bin path
        """
        raise NotImplementedError( "Should have implemented this method: '%s'" % self._assembler_bin.__name__ )


    def build_command_line(self, forward_fastq, reverse_fastq, workdir, read_correction='no', cpu=1, *args, **kwargs):
        """
        Build the command line to run (*args & **kwargs are potential parameters)
        The command must run in self.workdir with a number of cpu limited to self.cpu.
        It also has to take separate forward/reverse FastQ files as input.
        self.read_correction is a switch to take into account
        !!!! Important !!!!
        This method has to set the self.fasta_file which is the results of the assembly
        """

        raise NotImplementedError( "Should have implemented this method: '%s'" % self.build_command_line.__name__ )

    def run(self):
        """
        Check I/O and simply run the command
        """

        if self.cmd_line is None:
            logger.fatal('You have to call "build_command_line" method before "run" method')
            sys.exit("No command line available")

        for reads_file in (self.forward_fastq, self.reverse_fastq):
            if reads_file is None or not os.path.isfile(reads_file):
                logger.fatal('The input reads file does not exists:%s' % reads_file)
                sys.exit("Can't assemble %s" % reads_file)

        if self.workdir is None:
            logger.debug("Workdir is not set")
            sys.exit("Can't assemble %s" % self.forward_fastq)

        elif os.path.isdir(self.workdir):
            logger.debug("Remove previous assembler working dir before assembling:%s" % self.workdir)
            shutil.rmtree(self.workdir)
        os.mkdir(self.workdir)

        runner.logged_check_call(self.cmd_line)

        return self.fasta_file


class SGA(Assembler):

    def _assembler_wrapper(self):
        return binary_utils.Binary.assert_which('sga_assemble.py')

    def _assembler_bin(self):
        return binary_utils.Binary.assert_which('sga')

    def build_command_line(self, forward_fastq, reverse_fastq, workdir, read_correction='no', cpu=1, *args, **kwargs):
        self.forward_fastq = forward_fastq
        self.reverse_fastq = reverse_fastq
        self.workdir = workdir
        self.cpu = cpu
        self.read_correction = read_correction

        self.fasta_file =  os.path.join(workdir, 'assembly.fasta')
        logfile = os.path.join(workdir, 'assembly.log')
        tmp_dir = os.path.join(workdir, 'tmp')


        cmd_line = 'echo "component #' + self.forward_fastq + '" >> '
        cmd_line += logfile + ' && '
        cmd_line += self.assembler_wrapper + ' --forward ' + self.forward_fastq
        cmd_line += ' --reverse ' + self.reverse_fastq
        cmd_line += ' -o ' + self.fasta_file + ' --sga_bin ' + self.assembler_bin
        if self.read_correction in ('no', 'auto'):
            cmd_line += ' --no_correction' # !!! desactivate all SGA error corrections and filters
        cmd_line += ' --cpu ' + str(self.cpu)
        cmd_line += ' --tmp_dir %s' % tmp_dir
        cmd_line += ' >> ' + logfile + ' 2>&1'

        self.cmd_line = cmd_line

class AssemblerFactory:
    ASSEMBLER_ENGINES = [SGA,]

    def get(self, name):
        for assembler in self.ASSEMBLER_ENGINES:
            if assembler.name() == name:
                return assembler() #instantiate
        raise KeyError('Not a valid assembler. Valid keys: %s' % [a.name() for a in self.ASSEMBLER_ENGINES])
