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
import multiprocessing
import pysam
from docopt import docopt
import logging
import sys
from bitwise_flags import flags
from extract_pairs import ExtractPairs
from fastq_writer import FastqWriter
from read_evaluator import evaluate_pair
# from eval_worker import EvaluatorWorker as EW

__author__ = 'A. Jason Grundstad'


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
        self.fastq_queue_size = 256
        self.num_threads = num_threads
        self.FW = FastqWriter(outfile_stub=outfile)
        init_msg = """Input params -
        bamfile: {}
        outfiles: {}
        mouse_vars: {}
        fastq_queue_size: {}
        num_threads: {}"""
        self.logger.info(
            init_msg.format(bamfile, outfile + '_(1/2).fq.gz', mouse_vars,
                            self.fastq_queue_size, num_threads
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

    # def start_workers(self):
    #     for x in range(self.num_threads):
    #         worker = EW(queue=self.queue, number=x,
    #                     outfile=self.outfile)
    #         worker.daemon = True
    #         self.logger.info("Initiating worker {}".format(x))
    #         worker.start()

    def worker(self):
        read1, read2 = self.queue.get()
        evaluate_pair(read1=read1, read2=read2)

    def extract_reads(self):
        """
        Generator -
        Traverse .bam file for primary alignment pairs, skipping over
        secondary alignments.  pass primary alignments to evaluator and
        fill the fastq buffers.  dump to fq.gz files
        :return:
        """
        self.logger.info("Extracting primary-alignment read pairs.")

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
            yield (read1, read2)

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
        self.logger.info("Submitted {} sequence pairs to queue".format(
            pair_count))

    def fastq_queue_listener(self, queue):
        while True:
            print "HIHI"
            pair = queue.get()
            if pair == 'kill':
                break
            self.FW.print_reads(read1=pair[0], read2=pair[1])
            # self.logger.info("fastq_queue_listener: {}".format(pair[0].query_name))
            # self.FW.print_reads(pair[0], pair[1])

    def evaluate_reads(self):
        """
        Start fastq writing queue to assure paired reads are output properly
        Start the process pool over the readpair generator
        :return:
        """

        # for i,pair in enumerate(ExtractPairs(self.bam)):
        #     if i > 10:
        #         break
        #     else:
        #         print "{} - {}".format(pair[0].query_name, pair[1].query_name)
        manager = multiprocessing.Manager()
        fastq_queue = manager.Queue(maxsize=self.fastq_queue_size)
        fastq_writer_proc = multiprocessing.Process(
            target=self.fastq_queue_listener,
            args=(fastq_queue,))
        fastq_writer_proc.start()

        pool = multiprocessing.Pool(self.num_threads)

        # print "Popping 'hi' into the queue"
        # pool.apply_async(f, args=(fastq_queue, 'hi'))
        for i, pair in enumerate(ExtractPairs(self.bam)):
            if i % 100 == 0:
                self.logger.info("extracted {} pairs".format(i))
            pool.apply_async(evaluate_pair, args=(pair[0], pair[1],
                                                  fastq_queue))
        pool.close()
        pool.join()
        fastq_queue.put('kill')
        fastq_writer_proc.join()

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

def f(q, x):
    print "putting x: {}".format(x)
    q.put(x)

def main():
    args = docopt(__doc__)
    print args
    parser = BamParser(bamfile=args['-b'], outfile=args['-o'])
    parser.extract_reads()


if __name__ == '__main__':
    main()
