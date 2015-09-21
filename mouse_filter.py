"""
Remove all mouse-derived reads from a Human PDX .bam

Usage:
    mouse_filter.py -b BAMFILE -o OUT_FASTQ_STUB [-m MOUSE_VARS]

Options:
    -b BAMFILE           Human .bam alignment file
    -o OUT_FASTQ_STUB    Human fastq output file stub (e.g. 2015-123_human)
    -m MOUSE_VARS        Strain-specific variants
"""
import gzip
import multiprocessing
import pysam
from docopt import docopt
import logging
import sys
from bitwise_flags import flags
from extract_pairs import ExtractPairs
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
        self.fastq_queue_size = 100
        self.num_threads = num_threads
        self.outfile_stub = outfile
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

    def evaluate_reads(self):
        """
        Start fastq writing queue to assure paired reads are output properly
        Start the process pool over the readpair generator
        :return:
        """
        manager = multiprocessing.Manager()
        fastq_queue = manager.Queue(maxsize=self.fastq_queue_size)
        fastq_writer_proc = multiprocessing.Process(
            target=print_to_gzip,
            args=(fastq_queue, self.outfile_stub))
        fastq_writer_proc.start()

        pool = multiprocessing.Pool(self.num_threads)

        # pair = ExtractPairs(self.bam).next()
        # pool.apply_async(evaluate_pair, args=(read_to_dict(pair[0]),
        #                                       read_to_dict(pair[1]),
        #                                       fastq_queue))

        for i, pair in enumerate(ExtractPairs(self.bam)):
            if i % 1000000 == 0 and i > 0:
                # print pair
                self.logger.info("extracted {} pairs".format(1000000))
            pool.apply_async(evaluate_pair, args=(read_to_dict(pair[0]),
                                                  read_to_dict(pair[1]),
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

def read_to_dict(read):
    return {'query_name': read.query_name, 'seq': read.seq, 'qual': read.qual,
            'cigar': read.cigar, 'cigarstring': read.cigarstring,
            'flag': read.flag}

def read_to_fastq(read):
    return "@{}\n{}\n+\n{}\n".format(read.query_name, read.seq, read.qual)


def evaluate_pair(read1_dict, read2_dict, fastq_queue):
    """
    pop reads into queue if the reads are to be kept, pass if they are to be
    ignored.
    :param read1:
    :param read2:
    :return:
    """
    # print "putting dicts into queue: {}".format(read1_dict)
    # fastq_queue.put((read1_dict, read2_dict))
    if read1_dict['cigarstring'] and read2_dict['cigarstring']:
        # perfect alignment - pass over!
        if ((len(read1_dict['cigar']) == 1 and read1_dict['cigar'][0][0] ==0) and
                (len(read2_dict['cigar']) == 1 and read2_dict['cigar'][0][0] == 0)):
            pass
        else:
            fastq_queue.put((read1_dict, read2_dict))
    else:
        fastq_queue.put((read1_dict, read2_dict))

def read_to_fastq(read):
    return "@{}\n{}\n+\n{}\n".format(read['query_name'], read['seq'],
                                     read['qual'])

def print_to_gzip(queue, outfile_stub):
    fq1 = gzip.open(outfile_stub + '_1.fq.gz', 'w')
    fq2 = gzip.open(outfile_stub + '_2.fq.gz', 'w')
    while True:
        pair = queue.get()
        if pair == 'kill':
            print "returning 'kill'"
            break
        else:
            fq1.write(read_to_fastq(pair[0]))
            fq2.write(read_to_fastq(pair[1]))


def main():
    args = docopt(__doc__)
    print args
    parser = BamParser(bamfile=args['-b'], outfile=args['-o'])
    parser.evaluate_reads()


if __name__ == '__main__':
    main()
