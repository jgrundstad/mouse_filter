## .bam filtering tool for PDX samples
The goal is to remove all host DNA from a .bam, and regenerate the FASTQ
files for the standard alignment and variant calling pipeline. The assumed 
process leading up to this point:

1. Build custom model organism genome:
  1. Sequence normal tissue from the Xenograft model organism (e.g. Mouse).
  2. Generate a list of germline variants comparing the Mouse strain to mm10.
  3. Generate a custom reference by altering mm10 with the germline variants.
2. Align PDX sample reads to the custom reference
3. Isolate reads from all unaligned and imperfect alignments (using this tool!)
  * Output: 
    * Human_1.fq.gz
    * Human_2.fq.gz
    - ambiguous.bam
4. Investigate the ambiguous reads with Strain specific annotation, likely this will be added in to the Fastqs

```bash
usage: read_bam.py [-h] [-b BAM] -o OUTPUT [-c COMPRESSION]

Detect and isolate human reads from a bam file generated from human(SEQ)
aligned to mouse(REF). Accepts either: a file, or sam data piped from stdin.
NOTE: when reading from stdin, you must provide the SAM headers "@" via
samtools' -h flag.

optional arguments:
  -h, --help      show this help message and exit
  -b BAM          Input .bam (unsorted) [stdin]
  -o OUTPUT       Output stub e.g. Human.fastq
  -c COMPRESSION  Optional fq.gz compression rate [default: 4]
```

Note:  -m flag not implemented yet.