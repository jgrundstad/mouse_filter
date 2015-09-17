__author__ = 'jgrundst'


def evaluate_pair(read1=None, read2=None):
    """
    return True if the reads are to be kept, False if they are to be
    ignored.
    :param read1:
    :param read2:
    :return:
    """
    if read1.cigarstring and read2.cigarstring:
        if ((len(read1.cigar) == 1 and read1.cigar[0][0] ==0) and
                (len(read2.cigar) == 1 and read2.cigar[0][0] == 0)):
            pass
        else:
            return True
    else:
        return True
