"""
Remove all mouse-derived reads from a Human PDX .bam

Usage:
    mouse_filter.py -b BAMFILE -o OUT_FASTQ_STUB [-m MOUSE_VARS]

Options:
    -b BAMFILE           Human .bam alignment file
    -o OUT_FASTQ_STUB    Human fastq output file stub (e.g. 2015-123_human)
    -m MOUSE_VARS        Strain-specific variants
"""
from datetime import datetime
import pysam
from docopt import docopt
import logging
import sys
from bitwise_flags import flags
from fastq_writer import FastqWriter

__author__ = 'A. Jason Grundstad'


class BamParser:

    def __init__(self, bamfile=None, outfile=None, mouse_vars=None):
        self.handler = None
        self.logger = logging.getLogger(__name__)
        self.setup_logging()
        self.bam = None
        self.load_bam(bamfile)
        self.bamfilename = bamfile
        self.FW = FastqWriter(outfile_stub=outfile)
        init_msg = """Input params -
        bamfile: {}
        outfiles: {}
        mouse_vars: {}
        """
        self.logger.info(
            init_msg.format(bamfile, outfile + '_(1/2).fq.gz', mouse_vars))

    def setup_logging(self):
        if not self.logger.handlers:
            self.logger.setLevel(logging.DEBUG)
            self.handler = logging.FileHandler('MouseFilter.log')
            self.handler.setLevel(logging.DEBUG)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - [%(levelname)s]: %(message)s')
            self.handler.setFormatter(formatter)
            self.logger.addHandler(self.handler)
            self.logger.info("MouseFilter object created")

    def load_bam(self, bamfile):
        try:
            self.bam = pysam.AlignmentFile(bamfile, "rb")
        except IOError:
            self.logger.error("Cannot find bamfile: {}".format(bamfile))
            self.logger.info("Exiting")
            sys.exit(1)

    def extract_reads(self):
        """
        Traverse .bam file for primary alignment pairs, skipping over
        secondary alignments.  pass primary alignments to evaluator and
        fill the fastq buffers.  dump to fq.gz files
        :return:
        """
        self.logger.info("Extracting human reads to outfiles...")

        pair_count = 0
        keep_count = 0
        # the first read is always the primary alignment
        read1 = self.bam.next()
        read2 = None

        t1 = datetime.now()
        while read1:

            read2 = self.bam.next()
            # make sure we're looking at primary alignment of read2
            while ((bitwise_check(read2, 'second_in_pair') is False) and
                   (bitwise_check(read2, 'not_primary_alignment')) is True):
                read2 = self.bam.next()
            pair_count += 1
            if read1.query_name != read2.query_name:
                out_message = '''Pair #{}
                Read1 {}
                Read2 {}
                Don't have the same names.  Quitting...'''
                raise ValueError(out_message.format(pair_count, read1, read2))

            # eval the pair
            if evaluate_pair(read1=read1, read2=read2):
                keep_count += 1
                self.FW.print_reads(read1=read1, read2=read2)

            # Move on to next pair. Is there another read1?
            try:
                read1 = self.bam.next()
                while ((bitwise_check(read1, 'second_in_pair') is True)
                       and (bitwise_check(read1, 'not_primary_alignment')) is True):
                    read1 = self.bam.next()
            except StopIteration:
                self.logger.info('No more alignments, end of file{}'.format(
                    self.bamfilename
                ))
                read1 = None
        t2 = datetime.now()
        td = t2 - t1
        self.logger.info('Kept {} pairs out of {} in {}'.format(keep_count,
                                                                pair_count,
                                                                str(td)))

    def count_perfect_matches(self):
        c = 0
        for read in self.bam:
            if len(read.cigar) == 1 and read.cigar[0][0] == 0:
                c += 1
        print "{} perfectly aligned reads".format(c)


def bitwise_check(read, flag_string):
    if read.flag & flags[flag_string] == flags[flag_string]:
        return True
    else:
        return False

def read_to_fastq(read):
    return "@{}\n{}\n+\n{}\n".format(read.query_name, read.seq, read.qual)

def evaluate_pair(read1=None, read2=None):
    """
    return True if the reads are to be kept, False if they are to be
    ignored.
    :param read1:
    :param read2:
    :return:
    """
    if read1.cigarstring and read2.cigarstring:
        # perfect alignment - pass over!
        if ((len(read1.cigar) == 1 and read1.cigar[0][0] ==0) and
                (len(read2.cigar) == 1 and read2.cigar[0][0] == 0)):
            pass
        else:
            return True
    else:
        return True


def main():
    args = docopt(__doc__)
    parser = BamParser(bamfile=args['-b'], outfile=args['-o'])
    parser.extract_reads()

if __name__ == '__main__':
    main()
