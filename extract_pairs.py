import logging
from bitwise_flags import flags
__author__ = 'A. Jason Grundstad'

class ExtractPairs(object):
    """
    Generator to return pairs of primary read1/2 alignments from a
    pysam.AlignmentFile
    """

    def __init__(self, bam):
        """

        :param bam: pysam.AlignmentFile
        :return:
        """
        self.bam = bam
        self.pair_count = 0
        self.logger = logging.getLogger(__name__)
        self.read1 = self.bam.next()
        self.read2 = None

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        if self.read1:
            read2 = self.bam.next()
            while ((bitwise_flag_check(read2, 'second_in_pair') is False) and
                   (bitwise_flag_check(read2, 'not_primary_alignment')) is True):
                read2 = self.bam.next()
            self.read2 = read2
            self.pair_count += 1
            if self.read1.query_name != self.read2.query_name:
                out_message = '''Pair #{}
                Read1 {}
                Read2 {}
                Don't have the same names.  Quitting...'''
                raise ValueError(out_message.format(self.pair_count,
                                                    self.read1, self.read2))
            current_read1 = self.read1
            current_read2 = self.read2

            try:
                self.read1 = self.bam.next()
                while ((bitwise_flag_check(self.read1, 'second_in_pair') is True) and
                       (bitwise_flag_check(self.read1, 'not_primary_alignment')) is True):
                    self.read1 = self.bam.next()
            except StopIteration:
                self.logger.info('No more alignments, end of file{}'.format(
                    self.bam
                ))
                self.read1 = None
            return current_read1, current_read2
        else:
            raise StopIteration()


def bitwise_flag_check(read, flag_string):
    if read.flag & flags[flag_string] == flags[flag_string]:
        return True
    else:
        return False