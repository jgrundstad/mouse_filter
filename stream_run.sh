{ python read_bam.py -b test_data/1M.bam | pigz -p 3 -c - > Human_stream_1.fq.gz; } 2>&1 | pigz -p 3 -c - > Human_stream_2.fq.gz
