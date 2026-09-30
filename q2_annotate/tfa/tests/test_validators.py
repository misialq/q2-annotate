# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import re

import pandas as pd
import qiime2
from qiime2.plugin import Properties, ValidationError
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import (
    FeatureData,
    Taxonomy,
    TSVTaxonomyDirectoryFormat,
)

from q2_annotate.tfa.types._validators import validate_tfa_taxonomy


class TestTFATaxonomyValidator(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def setUp(self):
        super().setUp()
        self.taxonomy = pd.read_csv(
            self.get_data_path("gene-taxonomy.tsv"),
            sep="\t",
            index_col=0,
            dtype=str,
            keep_default_na=False,
        )
        self.large_taxonomy = pd.read_csv(
            self.get_data_path("gene-taxonomy-validation-levels.tsv"),
            sep="\t",
            index_col=0,
            dtype=str,
            keep_default_na=False,
        )

    def test_valid_mapping(self):
        """Accept a complete taxon/gene mapping at both validation levels."""
        for level in ("min", "max"):
            validate_tfa_taxonomy(self.taxonomy, level)

    def test_rejects_invalid_mapping(self):
        """Assert each TFA semantic validation message using invalid fixtures."""
        cases = {
            "gene-taxonomy-duplicate-id.tsv": (
                "TFA taxonomy feature IDs and taxon/gene pairs must be unique."
            ),
            "gene-taxonomy-duplicate-pair.tsv": (
                "TFA taxonomy feature IDs and taxon/gene pairs must be unique."
            ),
            "gene-taxonomy-empty-id.tsv": (
                "TFA taxonomy feature IDs, taxon IDs, and gene IDs must be nonempty."
            ),
            "gene-taxonomy-empty-taxon.tsv": (
                "TFA taxonomy taxon labels must be nonempty."
            ),
            "gene-taxonomy-bad-header.tsv": (
                "TFA taxonomy requires Taxon, Taxon ID, and Gene ID columns."
            ),
        }
        for filename, message in cases.items():
            invalid = pd.read_csv(
                self.get_data_path(filename),
                sep="\t",
                index_col=0,
                dtype=str,
                keep_default_na=False,
            )
            for level in ("min", "max"):
                with self.subTest(filename=filename, level=level):
                    with self.assertRaisesRegex(ValidationError, re.escape(message)):
                        validate_tfa_taxonomy(invalid, level)

    def test_uses_standard_taxonomy_format_and_transformers(self):
        """Round-trip all mapping columns through the standard taxonomy format."""
        artifact = qiime2.Artifact.import_data(
            "FeatureData[Taxonomy % Properties('tfa')]", self.taxonomy
        )
        self.assertEqual(artifact.format, TSVTaxonomyDirectoryFormat)
        self.assertEqual(artifact.type, FeatureData[Taxonomy % Properties("tfa")])
        self.assertLessEqual(artifact.type, FeatureData[Taxonomy])
        pd.testing.assert_frame_equal(artifact.view(pd.DataFrame), self.taxonomy)
        artifact.validate("max")

    def test_min_validation_ignores_errors_after_first_100_rows(self):
        """Check every row constraint beyond the minimal validation prefix."""
        cases = [
            (
                "Taxon ID",
                "",
                "TFA taxonomy feature IDs, taxon IDs, and gene IDs must be nonempty.",
            ),
            ("Taxon", "", "TFA taxonomy taxon labels must be nonempty."),
            (
                "Gene ID",
                self.large_taxonomy.iloc[0]["Gene ID"],
                "TFA taxonomy feature IDs and taxon/gene pairs must be unique.",
            ),
            (
                "Feature ID",
                self.large_taxonomy.index[0],
                "TFA taxonomy feature IDs and taxon/gene pairs must be unique.",
            ),
        ]
        for column, value, message in cases:
            with self.subTest(column=column):
                invalid = self.large_taxonomy.copy()
                if column == "Feature ID":
                    invalid = invalid.rename(index={invalid.index[-1]: value})
                else:
                    invalid.loc[invalid.index[-1], column] = value
                validate_tfa_taxonomy(invalid, "min")
                with self.assertRaisesRegex(ValidationError, re.escape(message)):
                    validate_tfa_taxonomy(invalid, "max")

    def test_extra_columns_are_preserved(self):
        """Allow additional annotations supported by the taxonomy format."""
        expected = self.taxonomy.assign(Source="annotation")
        artifact = qiime2.Artifact.import_data(
            "FeatureData[Taxonomy % Properties('tfa')]", expected
        )
        pd.testing.assert_frame_equal(artifact.view(pd.DataFrame), expected)

    def test_tfa_property_requires_pair_columns(self):
        """Require mapping columns only when the taxonomy has the tfa property."""
        ordinary = self.taxonomy[["Taxon"]]
        qiime2.Artifact.import_data("FeatureData[Taxonomy]", ordinary)
        with self.assertRaisesRegex(
            ValidationError,
            re.escape("TFA taxonomy requires Taxon, Taxon ID, and Gene ID columns."),
        ):
            qiime2.Artifact.import_data(
                "FeatureData[Taxonomy % Properties('tfa')]", ordinary
            )
