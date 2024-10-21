# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
import unittest
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt
from q2_moshpit.fastp import run_fastp

class TestFastp(unittest.TestCase):

    def setUp(self):
        self.input_sequences = SingleLanePerSampleSingleEndFastqDirFmt()
        self.output_sequences = SingleLanePerSampleSingleEndFastqDirFmt()

        # Create dummy input files
        for i in range(1, 4):
            input_fp = os.path.join(self.input_sequences.path, f'sample{i}.fastq.gz')
            with open(input_fp, 'w') as f:
                f.write('@SEQ_ID\nGATTTGGGGTTTCCCAGTTG\n+\nIIIIIIIIIIIIIIIIIIII\n')

    def test_run_fastp(self):
        output_sequences = run_fastp(self.input_sequences, trim_front1=5, trim_tail1=5, cut_window_size=4, cut_mean_quality=20, n_base_limit=5, length_required=15, qualified_quality_phred=15, unqualified_percent_limit=40.0, compression=2, thread=3)

        # Check if output files are created
        for i in range(1, 4):
            output_fp = os.path.join(output_sequences.path, f'sample{i}.fastq.gz')
            self.assertTrue(os.path.exists(output_fp))

        # Check if the output files are not empty
        for i in range(1, 4):
            output_fp = os.path.join(output_sequences.path, f'sample{i}.fastq.gz')
            self.assertGreater(os.path.getsize(output_fp), 0)

if __name__ == '__main__':
    unittest.main()
