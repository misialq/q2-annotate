# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import csv
import json
from pathlib import Path

import biom
import h5py
import numpy as np
import scipy.sparse as sp
from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.tfa import TFAFeatureTableFormat
from q2_annotate.tfa.types._transformer import (
    _biom_to_tfa_format,
    _tfa_format_to_biom,
)


class TestTFAFormat(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def _table_from_tsv(self):
        with open(
            self.get_data_path("tfa-table.tsv"), newline="", encoding="utf-8"
        ) as file:
            reader = csv.DictReader(file, delimiter="\t")
            sample_ids = reader.fieldnames[2:]
            rows = list(reader)

        return biom.Table(
            sp.csr_matrix(
                [[float(row[sample_id]) for sample_id in sample_ids] for row in rows]
            ),
            observation_ids=[
                json.dumps(
                    [row["taxon_id"], row["function_id"]],
                    ensure_ascii=False,
                    separators=(",", ":"),
                )
                for row in rows
            ],
            sample_ids=sample_ids,
        )

    def _write(self, path, table):
        with h5py.File(path, "w") as handle:
            table.to_hdf5(handle, generated_by="test")
        return TFAFeatureTableFormat(str(path), mode="r")

    def test_valid(self):
        path = Path(self.temp_dir.name) / "tfa-table.biom"
        self._write(path, self._table_from_tsv()).validate()

    def test_rejects_unpaired_id(self):
        table = self._table_from_tsv()
        observation_ids = list(table.ids(axis="observation"))
        observation_ids[0] = "T1"
        invalid = biom.Table(
            table.matrix_data.copy(),
            observation_ids=observation_ids,
            sample_ids=table.ids(axis="sample"),
        )
        path = Path(self.temp_dir.name) / "tfa-table.biom"
        with self.assertRaisesRegex(ValidationError, "pair ID"):
            self._write(path, invalid).validate()

    def test_rejects_negative_load(self):
        table = self._table_from_tsv()
        matrix = table.matrix_data.tolil()
        matrix[0, 0] = -1
        invalid = biom.Table(
            matrix.tocsr(),
            observation_ids=table.ids(axis="observation"),
            sample_ids=table.ids(axis="sample"),
        )
        path = Path(self.temp_dir.name) / "tfa-table.biom"
        with self.assertRaisesRegex(ValidationError, "nonnegative"):
            self._write(path, invalid).validate()

    def test_sparse_transformer_round_trip(self):
        expected = self._table_from_tsv()
        written = _biom_to_tfa_format(expected)
        readable = TFAFeatureTableFormat(str(written), mode="r")
        readable.validate()
        observed = _tfa_format_to_biom(readable)
        self.assertEqual(
            list(observed.ids(axis="sample")), list(expected.ids(axis="sample"))
        )
        self.assertEqual(
            list(observed.ids(axis="observation")),
            list(expected.ids(axis="observation")),
        )
        self.assertEqual(observed.matrix_data.nnz, 6)
        np.testing.assert_array_equal(
            observed.matrix_data.toarray(), expected.matrix_data.toarray()
        )
