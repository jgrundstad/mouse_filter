__author__ = 'jgrundst'

flags = {'read_paired': 1,
         'read_mapped_in_proper_pair': 2,
         'read_unmapped': 4,
         'mate_unmapped': 8,
         'read_reverse_strand': 16,
         'mate_reverse_strand': 32,
         'first_in_pair': 64,
         'second_in_pair': 128,
         'not_primary_alignment': 256,
         'read_fails_quality_checks': 512,
         'read_is_duplicate': 1024,
         'supplementary_alignment': 2048}
