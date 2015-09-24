import argparse
import gzip
import datetime
import pysam
__author__ = 'A. Jason Grundstad'

keep_count = 0
total_reads = 0

def read_bam(bamfile):
    global total_reads
    bam = pysam.AlignmentFile(bamfile, 'rb')
    read1 = bam.next()
    total_reads += 1
    read2 = None
    while read1:
        while read1.is_secondary and not read1.is_read1:
            read1 = bam.next()
            total_reads += 1
        read2 = bam.next()
        total_reads += 1
        while read2.is_secondary and not read2.is_read2:
            read2 = bam.next()
            total_reads += 1
        if read1.query_name is read2.query_name:
            print "{}\n{}".format(read1.query_name,
                                  read2.query_name)
            raise ValueError
        yield (read1, read2)

        read1 = bam.next()
        total_reads += 1
        while read1.is_secondary:
            read1 = bam.next()
            total_reads += 1


def print_fastq(outfile1=None, outfile2=None, read1=None, read2=None):
    outfile1.write("@{}\n{}\n+\n{}\n".format(read1.query_name, read1.seq,
                                             read1.qual))
    outfile2.write("@{}\n{}\n+\n{}\n".format(read2.query_name, read2.seq,
                                             read2.qual))
    global keep_count
    keep_count += 1


def main():
    desc='''
    Detect and isolate human reads from a bam file generated from human(SEQ)
    aligned to mouse(REF).  Accepts either: a file, or sam data piped from stdin.
    NOTE: when reading from stdin, you must provide the SAM headers "@" via
    samtools' -h flag.
    '''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-b', dest='bam', default='-',
                        help="Input .bam (unsorted) [stdin]")
    parser.add_argument('-o', dest='output', required=True,
                        help='Output stub e.g. Human.fastq')
    parser.add_argument('-c', dest='compression', required=False, default=4,
                        type=int,
                        help='Optional fq.gz compression rate [default: 4]')
    args = parser.parse_args()

    fq1 = gzip.open(args.output + '_1.fq.gz', 'w',
                    compresslevel=args.compression)
    fq2 = gzip.open(args.output + '_2.fq.gz', 'w',
                    compresslevel=args.compression)

    pair_count = None
    t_start = datetime.datetime.now()

    for pair_count, pair in enumerate(read_bam(args.bam), start=1):
        # do we have alignments
        if not (pair[0].is_unmapped and pair[1].is_unmapped):
            # are both alignments perfect, pass it over
            if ((len(pair[0].cigar) == 1 and pair[0].cigar[0][0] == 0) and
                    (len(pair[1].cigar) == 1 and pair[1].cigar[0][0] == 0)):
                pass
            else:
                print_fastq(outfile1=fq1, outfile2=fq2,
                            read1=pair[0], read2=pair[1])
        else:
            print_fastq(outfile1=fq1, outfile2=fq2,
                        read1=pair[0], read2=pair[1])

    t_end = datetime.datetime.now()
    time_delta = t_end - t_start
    pct = (keep_count + 0.0) / pair_count * 100
    print "{} total reads".format(total_reads)
    print "kept {} primary alignment pairs out of {}  %{:.4f}".format(
        keep_count, pair_count, pct)
    print "time: {}".format(str(time_delta))

if __name__ == '__main__':
    main()
