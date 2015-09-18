__author__ = 'A. Jason Grundstad'


def evaluate_pair(read1, read2, fastq_queue):
    """
    pop reads into queue if the reads are to be kept, pass if they are to be
    ignored.
    :param read1:
    :param read2:
    :return:
    """
    fastq_queue.put(('strang1', 'strang2'))
    return True
    print "{}".format(read1.query_name)
    if read1.cigarstring and read2.cigarstring:
        # perfect alignment - pass over!
        if ((len(read1.cigar) == 1 and read1.cigar[0][0] ==0) and
                (len(read2.cigar) == 1 and read2.cigar[0][0] == 0)):
            pass
        else:
            fastq_queue.put((read1, read2))
    else:
        fastq_queue.put((read1, read2))
