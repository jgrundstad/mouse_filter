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