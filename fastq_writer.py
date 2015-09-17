import gzip
import logging
import sys

__author__ = 'jgrundst'


class FastqWriter:

    def __init__(self, outfile_stub=None):
        self.logger = logging.getLogger(__name__)
        self.outfile_stub = outfile_stub
        self.fq1 = None
        self.fq2 = None
        self.set_outfiles()

    def set_outfiles(self):
        try:
            self.fq1 = gzip.open(self.outfile_stub + '_1.fq.gz', 'w')
            self.fq2 = gzip.open(self.outfile_stub + '_2.fq.gz', 'w')
        except IOError:
            self.logger.error(
                "Unable to open outfile for writing: {}".format(
                    self.outfile_stub
                )
            )
            self.logger.critical("Exiting")
            sys.exit(1)

    def print_reads(self, read1=None, read2=None):
        self.fq1.write(read_to_fastq(read1))
        self.fq2.write(read_to_fastq(read2))

def read_to_fastq(read):
    return "@{}\n{}\n+\n{}\n".format(read.query_name, read.seq, read.qual)
