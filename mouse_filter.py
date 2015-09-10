"""
Remove all mouse-derived reads from a Human PDX .bam

Usage:
    mouse_filter.py -b BAMFILE -o OUT_FASTQ [-m MOUSE_VARS]

Options:
    -b BAMFILE      Human .bam alignemnt file
    -o OUT_FASTQ    Human fastq output file
    -m MOUSE_VARS   Strain-specific variants
"""
__author__ = 'jgrundst'
import pysam
from docopt import docopt
import logging
import sys


class MouseFilter:

    def __init__(self, bamfile=None, outfile=None, mouse_vars=None):
        self.handler = None
        self.logger = logging.getLogger(__name__)
        self.setup_logging()
        self.bam = None
        self.load_bam(bamfile)
        self.outfile = None
        self.set_outfile(outfile)
        self.headers = None

        self.logger.info(
            "Input params - \nbamfile: {}\noutfile: {}\nmouse_vars: {}".format(
                bamfile, outfile, mouse_vars
            )
        )

    def load_bam(self, bamfile):
        try:
            self.bam = pysam.AlignmentFile(bamfile, "rb")
        except IOError:
            self.logger.error("Cannot find bamfile: {}".format(bamfile))
            self.logger.info("Exiting")
            sys.exit(1)

    def set_outfile(self, outfile):
        try:
            self.outfile = open(outfile, 'w')
        except IOError:
            self.logger.error(
                "Unable to open outfile for writing: {}".format(
                    outfile
                )
            )
            self.logger.info("Exiting")
            sys.exit(1)

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

    def show_headers(self):
        self.headers = self.bam.header
        for record_type, records in self.headers.items():
            print record_type
            for i, record in enumerate(records):
                print "\t{},".format(i+1)
                if type(record) == str:
                    print "\t\t{}".format(record)
                elif type(record) == dict:
                    for field, value in record.items():
                        print "\t\t{}\t{}".format(field, value)



def main():
    args = docopt(__doc__)
    print args
    mf = MouseFilter(bamfile=args['-b'], outfile=args['-o'])
    mf.show_headers()



if __name__ == '__main__':
    main()