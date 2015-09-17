"""
Remove all mouse-derived reads from a Human PDX .bam

Usage:
    mouse_filter.py -b BAMFILE -o OUT_FASTQ_STUB [-m MOUSE_VARS]

Options:
    -b BAMFILE           Human .bam alignment file
    -o OUT_FASTQ_STUB    Human fastq output file stub (e.g. 2015-123_human)
    -m MOUSE_VARS        Strain-specific variants
"""
import Queue
import pysam
from docopt import docopt
import logging
import sys
from bitwise_flags import flags
from eval_worker import EvaluatorWorker as EW
from fastq_writer import FastqWriter

__author__ = 'jgrundst'


class BamParser:

    def __init__(self, bamfile=None, outfile=None, mouse_vars=None,
                 num_threads=4, queue_maxsize=256):
        self.handler = None
        self.logger = logging.getLogger(__name__)
        self.setup_logging()
        self.bam = None
        self.load_bam(bamfile)
        self.bamfilename = bamfile
        self.headers = None
        self.queue = Queue.Queue(maxsize=queue_maxsize)
        self.FW = FastqWriter(outfile_stub=outfile)
        self.num_threads = num_threads
        init_msg = """Input params -
        bamfile: {}
        outfiles: {}
        mouse_vars: {}
        queue_maxsize: {}
        num_threads: {}"""
        self.logger.info(
            init_msg.format(bamfile, outfile + '_(1/2).fq.gz', mouse_vars,
                            queue_maxsize, num_threads
            )
        )

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

    def start_workers(self):
        for x in range(self.num_threads):
            worker = EW(FW=self.FW, queue=self.queue)
            worker.daemon = True
            self.logger.info("Initiating worker {}".format(x))
            worker.start()

    def extract_reads(self):
        """
        Traverse .bam file for primary alignment pairs, skipping over
        secondary alignments.  pass primary alignments to evaluator and
        fill the fastq buffers.  dump to fq.gz files
        :return:
        """
        self.logger.info("Extracting human reads to outfiles...")
        self.logger.info("Initializing queue with {} threads".format(
            self.num_threads
        ))

        self.start_workers()

        pair_count = 0
        # the first read is always the primary alignment
        read1 = self.bam.next()
        read2 = None
        while read1:

            read2 = self.bam.next()
            # make sure we're looking at primary alignment of read2
            while ((bitwise_flag_check(read2, 'second_in_pair') is False) and
                   (bitwise_flag_check(read2, 'not_primary_alignment')) is True):
                read2 = self.bam.next()
            pair_count += 1
            if read1.query_name != read2.query_name:
                out_message = '''Pair #{}
                Read1 {}
                Read2 {}
                Don't have the same names.  Quitting...'''
                raise ValueError(out_message.format(pair_count, read1, read2))

            # queue the evaluation job for the workers
            self.queue.put((read1, read2))

            # Move on to next pair. Is there another read1?
            try:
                read1 = self.bam.next()
                while ((bitwise_flag_check(read1, 'second_in_pair') is True) and
                       (bitwise_flag_check(read1, 'not_primary_alignment')) is True):
                    read1 = self.bam.next()
            except StopIteration:
                self.logger.info('No more alignments, end of file{}'.format(
                    self.bamfilename
                ))
                read1 = None
        # put final pair in queue
        self.queue.put((read1, read2))
        self.queue.join()
        self.logger.info("Submitted {} sequence pairs to queue".format(
            pair_count))

    def count_perfect_matches(self):
        c = 0
        for read in self.bam:
            if len(read.cigar) == 1 and read.cigar[0][0] == 0:
                c += 1
        print "{} perfectly aligned reads".format(c)


def bitwise_flag_check(read, flag_string):
    if read.flag & flags[flag_string] == flags[flag_string]:
        return True
    else:
        return False


def read_to_fastq(read):
    return "@{}\n{}\n+\n{}\n".format(read.query_name, read.seq, read.qual)


def main():
    args = docopt(__doc__)
    print args
    parser = BamParser(bamfile=args['-b'], outfile=args['-o'])
    parser.extract_reads()

if __name__ == '__main__':
    main()
