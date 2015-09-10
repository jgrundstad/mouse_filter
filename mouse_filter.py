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

class MouseFilter():

    def __init__(self, bamfile=None, outfile=None, mouse_vars=None):
        self.setup_logging()

        try:
            self.bamfile = pysam.AlignmentFile(bamfile, "rb")
        except IOError:
            self.logger.error("Cannot find bamfile: {}".format(bamfile))
            self.logger.info("Exiting")
            sys.exit(1)

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

        self.outfile = outfile
        self.logger.info(
            "Input params - \nbamfile: {}\noutfile: {}\nmouse_vars: {}".format(
                bamfile, outfile, mouse_vars
            )
        )

    def setup_logging(self):
        self.logger = logging.getLogger(__name__)
        if not self.logger.handlers:
            self.logger.setLevel(logging.DEBUG)
            self.handler = logging.FileHandler('MouseFilter.log')
            self.handler.setLevel(logging.DEBUG)
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - [%(levelname)s]: %(message)s')
            self.handler.setFormatter(formatter)
            self.logger.addHandler(self.handler)
            self.logger.info("MouseFilter object created")

    def extract_human_reads(self):



def main():
    args = docopt(__doc__)
    print args
    mf = MouseFilter(bamfile=args['-b'], outfile=args['-o'])
    mf.extract_human_reads()


if __name__ == '__main__':
    main()