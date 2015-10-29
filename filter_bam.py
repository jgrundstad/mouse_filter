import argparse
import gzip
import datetime
import pysam
import sys

__author__ = 'A. Jason Grundstad'

total_reads = 0
window = 500


def read_bam(bamfile):
    global total_reads
    bam = pysam.AlignmentFile(bamfile, 'rb')
    read1 = bam.next()
    total_reads += 1
    read2 = None
    while read1:
        while read1.is_secondary and not read1.is_read1 and not read1.is_supplementary:
            read1 = bam.next()
            total_reads += 1
        read2 = bam.next()
        total_reads += 1
        while read2.is_secondary and not read2.is_read2 and not read2.is_supplementary:
            read2 = bam.next()
            total_reads += 1
        if read1.query_name != read2.query_name:
            print "{}\n{}".format(read1.query_name,
                                  read2.query_name)
            raise ValueError
        yield (read1, read2)

        read1 = bam.next()
        total_reads += 1
        while read1.is_secondary:
            read1 = bam.next()
            total_reads += 1


def print_fastq_to_pipes(read1=None, read2=None, **kwargs):
    sys.stdout.write("@{}\n{}\n+\n{}\n".format(read1.query_name, read1.seq,
                                               read1.qual))
    sys.stderr.write("@{}\n{}\n+\n{}\n".format(read2.query_name, read2.seq,
                                               read2.qual))


def print_fastq(outfile1=None, outfile2=None, read1=None, read2=None):
    outfile1.write("@{}\n{}\n+\n{}\n".format(read1.query_name, read1.seq,
                                             read1.qual))
    outfile2.write("@{}\n{}\n+\n{}\n".format(read2.query_name, read2.seq,
                                             read2.qual))


def both_mapped(read1, read2):
    return not (read1.is_unmapped and read2.is_unmapped)


def mated(read1, read2):
    return (read1.rname == read2.rname) and (abs(read1.isize) < window)


def perfect_alignments(read1, read2):
    return ((len(read1.cigar) == 1 and read1.cigar[0][0] == 0) and
            (len(read2.cigar) == 1 and read2.cigar[0][0] == 0))


def evaluate(bam=None, fq1=None, fq2=None):
    pair_count = 0
    keep_count = 0
    ambiguous_count = 0
    print_it = print_fastq
    if fq1 is None:
        print_it = print_fastq_to_pipes

    for pair_count, pair in enumerate(read_bam(bam), start=1):
        # do we have alignments
        if not (pair[0].is_unmapped and pair[1].is_unmapped):
            ''' are both alignments perfect, pass it over unless:
            *  read1.reference_id != read2.reference_id
            *  NM in read.tags indicates edit distance from reference
            *  insert size: read1.isize , negative for read2
            '''
            if (both_mapped(pair[0], pair[1]) and
                    mated(pair[0], pair[1]) and
                    perfect_alignments(pair[0], pair[1])):
                pass
            else:
                #  Consider sending imperfect alignments to other pair of files
                print_it(outfile1=fq1, outfile2=fq2, read1=pair[0], read2=pair[1])
                ambiguous_count += 1
        else:
            print_it(outfile1=fq1, outfile2=fq2, read1=pair[0], read2=pair[1])
            keep_count += 1
    return pair_count, keep_count, ambiguous_count


def main():
    desc = '''
    Detect and isolate human reads from a bam file generated from human(SEQ)
    aligned to mouse(REF).  Accepts either: a file, or sam data piped from stdin.
    NOTE: when reading from stdin, you must provide the SAM headers "@" via
    samtools' -h flag.
    '''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-b', dest='bam', default='-',
                        help="Input .bam (unsorted) [stdin]")
    parser.add_argument('-o', dest='output', required=False,
                        help='Output stub e.g. Human.fastq')
    parser.add_argument('-c', dest='compression', required=False, default=4,
                        type=int,
                        help='Optional fq.gz compression rate [default: 4]')
    args = parser.parse_args()

    logfile = open('runlog.txt', 'a')
    fq1 = None
    fq2 = None
    if args.output:
        fq1 = gzip.open(args.output + '_1.fq.gz', 'w',
                        compresslevel=args.compression)
        fq2 = gzip.open(args.output + '_2.fq.gz', 'w',
                        compresslevel=args.compression)

    t_start = datetime.datetime.now()
    print >>logfile, "----------\nStart time: {}".format(t_start)
    pair_count, keep_count, ambiguous_count = evaluate(
        fq1=fq1, fq2=fq2, bam=args.bam)
    t_end = datetime.datetime.now()
    print >>logfile, "End time:   {}".format(t_end)

    time_delta = t_end - t_start
    keep_pct = (keep_count + 0.0) / pair_count * 100
    ambig_pct = (ambiguous_count + 0.0) / pair_count * 100
    global total_reads
    print >>logfile, "{} total reads".format(total_reads)
    print >>logfile, "kept {} alignment pairs out of {}  %{:.4f}".format(
        keep_count, pair_count, keep_pct)
    print >>logfile, "kept {} ambiguous alignment pairs out of {}  %{:.4f}".format(
        ambiguous_count, pair_count, ambig_pct)
    print >>logfile, "time delta: {}".format(str(time_delta))

if __name__ == '__main__':
    main()
