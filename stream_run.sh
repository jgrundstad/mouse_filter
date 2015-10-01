#!/usr/bin/env bash
function helptext {
    echo "-b    Bam file to be processed"
    echo "-o    output file stub"
    echo "-h    this message"
    exit 0
}
[[ $# -gt 0 ]] || { helptext; }

while getopts o:b:h OPT; do
    case $OPT in
        o)
            OUTSTUB=$OPT
            ;;
        b)
            BAM=$OPT
            ;;
        h)
            helptext
            ;;
        \?)
            echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
            helptext
            ;;
    esac
done

{ python read_bam.py -b $BAM | gzip -4 -c - > $OUTSTUB_1.fq.gz; } 2>&1 | gzip -4 -c - > $OUTSTUB_2.fq.gz
