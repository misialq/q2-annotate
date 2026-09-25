# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import tempfile
import unittest
from pathlib import Path

import biom
import h5py
import numpy as np
from qiime2.plugin import ValidationError

from q2_annotate.tfa import TFAFeatureTableFormat
from q2_annotate.tfa.types._transformer import (
    _biom_to_tfa_format,
    _tfa_format_to_biom,
)


class TestTFAFormat(unittest.TestCase):
    def _write(self, path, observation_id, value):
        table = biom.Table(
            np.array([[value]]), observation_ids=[observation_id], sample_ids=["S1"]
        )
        with h5py.File(path, "w") as handle:
            table.to_hdf5(handle, generated_by="test")
        return TFAFeatureTableFormat(path, mode="r")

    def test_valid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "tfa-table.biom"
            self._write(path, json.dumps(["T1", "gene"]), 1.25).validate()

    def test_rejects_unpaired_id(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "tfa-table.biom"
            with self.assertRaisesRegex(ValidationError, "pair ID"):
                self._write(path, "T1", 1).validate()

    def test_rejects_negative_load(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "tfa-table.biom"
            with self.assertRaisesRegex(ValidationError, "nonnegative"):
                self._write(path, json.dumps(["T1", "gene"]), -1).validate()

    def test_sparse_transformer_round_trip(self):
        pair_id = json.dumps(["T1", "gene"])
        expected = biom.Table(
            np.array([[1.25, 0, 0]]),
            observation_ids=[pair_id],
            sample_ids=["S1", "S2", "S3"],
        )
        written = _biom_to_tfa_format(expected)
        readable = TFAFeatureTableFormat(str(written), mode="r")
        readable.validate()
        observed = _tfa_format_to_biom(readable)
        self.assertEqual(list(observed.ids(axis="sample")), ["S1", "S2", "S3"])
        self.assertEqual(list(observed.ids(axis="observation")), [pair_id])
        self.assertEqual(observed.matrix_data.nnz, 1)
        np.testing.assert_array_equal(observed.matrix_data.toarray(), [[1.25, 0, 0]])
