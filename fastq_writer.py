import gzip
import logging
import sys

__author__ = 'jgrundst'


class FastqWriter:

    def __init__(self, outfile_stub=None, paired=None):
        self.logger = logging.getLogger(__name__)
        self.outfile_stub = outfile_stub
        self.paired = paired
        self.fq1 = None
        self.fq1_buffer = None
        if paired is True:
            self.fq2 = None
            self.fq2_buffer = None
        self.buffer_read_count = None
        self.buffer_size = 10

        self.set_outfiles()

    def set_outfiles(self):
        try:
            self.fq1 = gzip.open(self.outfile_stub + '_1.fq.gz', 'w')
            if self.paired:
                self.fq2 = gzip.open(self.outfile_stub + '_2.fq.gz', 'w')
        except IOError:
            self.logger.error(
                "Unable to open outfile for writing: {}".format(
                    self.outfile_stub
                )
            )
            self.logger.critical("Exiting")
            sys.exit(1)

    def buffer_reads(self, read1=None, read2=None):
        self.fq1_buffer += read_to_fastq(read1)
        if read2:
            self.fq2_buffer += read_to_fastq(read2)
        self.buffer_read_count += 1
        if self.buffer_read_count == self.buffer_size:
            self.print_fq_buffers()

    def print_fq_buffers(self):
        self.fq1.write(self.fq1_buffer)
        self.fq1_buffer = ''
        if self.fq2:
            self.fq2.write(self.fq2_buffer)
            self.fq2_buffer = ''
        self.buffer_read_count = 0

    def flush_buffers(self):
        self.print_fq_buffers()
        self.fq1.close()
        self.fq2.close()


def read_to_fastq(read):
    return "@{}\n{}\n+\n{}\n".format(read.query_name, read.seq, read.qual)
