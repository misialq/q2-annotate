# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import biom
import numpy as np
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.plugin_setup import plugin


class TestTFAAction(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def test_action_registration(self):
        self.assertIn("estimate_tfa", plugin.methods)
        self.assertNotIn("_estimate_tfa_table", plugin.methods)
        self.assertNotIn("tfa_by_taxon", plugin.methods)
        self.assertNotIn("tfa_by_function", plugin.methods)
        self.assertNotIn("estimate_tfa", plugin.pipelines)
        self.assertNotIn("TFA", plugin.type_fragments)

    def _inputs(self, inventory_contigs=("C1", "C2", "C3")):
        abundance = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]",
            biom.Table(
                np.array([[2, 0], [1, 3], [0, 4]]),
                observation_ids=["C1", "C2", "C3"],
                sample_ids=["S1", "S2"],
            ),
        )
        inventory = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]",
            biom.Table(
                np.array([[2, 1, 0], [0, 1, 3]])[:, : len(inventory_contigs)],
                observation_ids=["geneA", "gene B"],
                sample_ids=list(inventory_contigs),
            ),
        )
        taxonomy = qiime2.Artifact.import_data(
            "FeatureData[Taxonomy]",
            pd.DataFrame(
                {"Taxon": ["k__Bacteria; p__A", "k__Bacteria; p__B"]},
                index=pd.Index(["T1", "T2"], name="Feature ID"),
            ),
        )
        mapping = qiime2.Artifact.import_data(
            "FeatureMap[TaxonomyToContigs]",
            {"T1": ["C1", "C2"], "T2": ["C3"]},
        )
        return dict(
            abundance_matrix=abundance,
            feature_inventory=inventory,
            taxonomy=taxonomy,
            taxon_to_contig_map=mapping,
        )

    def test_estimate_returns_frequency_and_gene_taxonomy(self):
        result = plugin.methods["estimate_tfa"](**self._inputs())
        self.assertEqual(len(result), 2)
        self.assertEqual(str(result.feature_load.type), "FeatureTable[Frequency]")
        self.assertEqual(str(result.gene_taxonomy.type), "FeatureData[GeneTaxonomy]")
        self.assertEqual(
            list(result.feature_load.view(biom.Table).ids(axis="sample")),
            ["S1", "S2"],
        )
        self.assertEqual(
            list(result.feature_load.view(biom.Table).ids(axis="observation")),
            list(result.gene_taxonomy.view(pd.DataFrame).index),
        )

    def test_empty_table_when_no_contigs_overlap(self):
        result = plugin.methods["estimate_tfa"](**self._inputs(("C4",)))
        self.assertEqual(result.feature_load.view(biom.Table).shape, (0, 2))
        self.assertEqual(result.gene_taxonomy.view(pd.DataFrame).shape, (0, 3))
