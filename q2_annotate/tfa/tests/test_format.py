# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
import qiime2
from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.tfa import GeneTaxonomyFormat
from q2_annotate.tfa.types._transformer import (
    _gene_taxonomy_to_format,
    _gene_taxonomy_to_dataframe,
)


class TestGeneTaxonomyFormat(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def test_valid_and_empty_mapping(self):
        for filename in ("gene-taxonomy.tsv", "gene-taxonomy-empty.tsv"):
            GeneTaxonomyFormat(self.get_data_path(filename), mode="r").validate("max")

    def test_rejects_invalid_mapping(self):
        for filename in (
            "gene-taxonomy-duplicate-id.tsv",
            "gene-taxonomy-duplicate-pair.tsv",
            "gene-taxonomy-empty-id.tsv",
            "gene-taxonomy-bad-header.tsv",
        ):
            with self.assertRaises(ValidationError):
                GeneTaxonomyFormat(self.get_data_path(filename), mode="r").validate(
                    "max"
                )

    def test_mapping_transformer_round_trip_preserves_ids_and_empty_labels(self):
        expected = pd.read_csv(
            self.get_data_path("gene-taxonomy.tsv"),
            sep="\t",
            index_col=0,
            dtype=str,
            keep_default_na=False,
        )
        written = _gene_taxonomy_to_format(expected)
        readable = GeneTaxonomyFormat(str(written), mode="r")
        readable.validate("max")
        pd.testing.assert_frame_equal(_gene_taxonomy_to_dataframe(readable), expected)
        artifact = qiime2.Artifact.import_data("FeatureData[GeneTaxonomy]", expected)
        pd.testing.assert_frame_equal(artifact.view(pd.DataFrame), expected)

    def test_transformer_rejects_duplicate_ids(self):
        invalid = pd.read_csv(
            self.get_data_path("gene-taxonomy-duplicate-id.tsv"),
            sep="\t",
            index_col=0,
            dtype=str,
            keep_default_na=False,
        )
        with self.assertRaisesRegex(ValueError, "unique"):
            _gene_taxonomy_to_format(invalid)
