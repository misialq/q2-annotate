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
import scipy.sparse as sp
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.tfa import estimate_tfa
from q2_annotate.tfa.types._utils import _gene_taxonomy_id


class TestTFA(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def setUp(self):
        super().setUp()
        self.abundance = pd.read_csv(
            self.get_data_path("contig-abundance.tsv"), sep="\t", index_col=0
        ).T
        self.inventory = pd.read_csv(
            self.get_data_path("gene-inventory.tsv"), sep="\t", index_col=0
        ).T
        self.taxonomy = pd.read_csv(
            self.get_data_path("estimate-taxonomy.tsv"), sep="\t", index_col=0
        )
        mapping = pd.read_csv(self.get_data_path("taxon-to-contigs.tsv"), sep="\t")
        self.mapping = mapping.groupby("Taxon ID")["Contig ID"].apply(list).to_dict()

    def estimate(self, abundance=None, inventory=None, taxonomy=None, mapping=None):
        abundance = self.abundance if abundance is None else abundance
        inventory = self.inventory if inventory is None else inventory
        taxonomy = self.taxonomy if taxonomy is None else taxonomy
        mapping = self.mapping if mapping is None else mapping
        abundance, inventory = abundance.T, inventory.T
        return estimate_tfa(
            biom.Table(
                abundance.values,
                observation_ids=abundance.index,
                sample_ids=abundance.columns,
            ),
            biom.Table(
                inventory.values,
                observation_ids=inventory.index,
                sample_ids=inventory.columns,
            ),
            taxonomy,
            mapping,
        )

    def test_sample_resolved_values_and_sparse_pairs(self):
        observed, gene_taxonomy = self.estimate()
        self.assertEqual(list(observed.ids(axis="sample")), ["S1", "S2", "S3"])
        self.assertEqual(observed.shape, (4, 3))
        self.assertTrue(sp.issparse(observed.matrix_data))
        self.assertEqual(observed.matrix_data.nnz, 12)
        expected = {
            ("T1", "bla_TEM"): [200, 240, 180],
            ("T1", "vanA"): [100, 120, 90],
            ("T2", "bla_TEM"): [50, 40, 60],
            ("T2", "mecA"): [150, 120, 180],
        }
        for pair_id, values in zip(
            observed.ids(axis="observation"), observed.matrix_data.toarray()
        ):
            np.testing.assert_array_equal(
                values,
                expected[tuple(gene_taxonomy.loc[pair_id, ["Taxon ID", "Gene ID"]])],
            )

    def test_no_common_contigs_retains_samples(self):
        observed, gene_taxonomy = self.estimate(mapping={"T1": ["not-a-contig"]})
        self.assertEqual(observed.shape, (0, 3))
        self.assertEqual(list(observed.ids(axis="sample")), ["S1", "S2", "S3"])

    def test_ignores_extra_contigs(self):
        abundance = self.abundance.assign(C7=[1000, 1000, 1000])
        taxonomy = pd.DataFrame(
            {"Taxon": ["E. coli", "K. pneumoniae", "Other"]},
            index=["T1", "T2", "T3"],
        )
        mapping = {**self.mapping, "T3": ["C7", "C8"]}
        observed, gene_taxonomy = self.estimate(
            abundance=abundance, taxonomy=taxonomy, mapping=mapping
        )
        self.assertEqual(observed.shape, (4, 3))
        self.assertEqual(observed.matrix_data.nnz, 12)

    def test_zero_abundance_keeps_observed_pairs(self):
        observed, gene_taxonomy = self.estimate(abundance=self.abundance * 0)
        self.assertEqual(observed.shape, (4, 3))
        self.assertEqual(observed.matrix_data.nnz, 0)

    def test_extreme_abundances(self):
        abundance = self.abundance * 0.0
        abundance.loc["S1", "C1"] = 1e-9
        abundance.loc["S1", "C2"] = 1e12
        observed, gene_taxonomy = self.estimate(abundance=abundance)
        loads = {
            tuple(gene_taxonomy.loc[pair_id, ["Taxon ID", "Gene ID"]]): values
            for pair_id, values in zip(
                observed.ids(axis="observation"), observed.matrix_data.toarray()
            )
        }
        np.testing.assert_allclose(
            loads[("T1", "bla_TEM")],
            [2e-9, 0, 0],
            rtol=1e-12,
            atol=0,
        )
        np.testing.assert_allclose(
            loads[("T2", "mecA")],
            [3e12, 0, 0],
            rtol=1e-12,
            atol=0,
        )

    def test_mapping_preserves_ambiguous_characters(self):
        inventory = self.inventory.rename(columns={"bla_TEM": 'a,"b'})
        taxonomy = self.taxonomy.rename(index={"T1": 't,"1'})
        mapping = {'t,"1': self.mapping["T1"], "T2": self.mapping["T2"]}
        observed, gene_taxonomy = self.estimate(
            inventory=inventory, taxonomy=taxonomy, mapping=mapping
        )
        pairs = set(
            gene_taxonomy[["Taxon ID", "Gene ID"]].itertuples(index=False, name=None)
        )
        self.assertIn(('t,"1', 'a,"b'), pairs)

    def test_ambiguous_taxon_assignment_rejected(self):
        mapping = {**self.mapping, "T2": ["C1", "C2", "C4", "C6"]}
        with self.assertRaisesRegex(ValueError, "multiple taxa"):
            self.estimate(mapping=mapping)

    def test_stable_feature_ids_and_mapping(self):
        table, mapping = self.estimate()
        reordered, reordered_mapping = self.estimate(
            inventory=self.inventory[self.inventory.columns[::-1]],
            abundance=self.abundance[self.abundance.columns[::-1]],
        )
        pd.testing.assert_frame_equal(
            mapping.sort_index(), reordered_mapping.sort_index()
        )
        self.assertEqual(list(table.ids(axis="observation")), list(mapping.index))
        self.assertEqual(
            set(mapping.index),
            {"T1|bla_TEM", "T1|vanA", "T2|bla_TEM", "T2|mecA"},
        )
        pd.testing.assert_frame_equal(
            table.to_dataframe(dense=True).sort_index(),
            reordered.to_dataframe(dense=True).sort_index(),
        )

    def test_readable_ids_escape_ambiguous_components(self):
        cases = pd.read_csv(self.get_data_path("pair-ids.tsv"), sep="\t", dtype=str)
        observed = [
            _gene_taxonomy_id(taxon, gene)
            for taxon, gene in zip(cases["Taxon ID"], cases["Gene ID"])
        ]
        self.assertEqual(observed, list(cases["Feature ID"]))
        self.assertEqual(len(observed), len(set(observed)))

    def test_registered_estimate_outputs_and_fractional_values(self):
        import qiime2
        import biom

        abundance = self.abundance.T * 0.01
        inventory = self.inventory.T
        result = self.plugin.methods["estimate_tfa"](
            abundance_matrix=qiime2.Artifact.import_data(
                "FeatureTable[Frequency]",
                biom.Table(
                    abundance.values,
                    observation_ids=abundance.index,
                    sample_ids=abundance.columns,
                ),
            ),
            feature_inventory=qiime2.Artifact.import_data(
                "FeatureTable[Frequency]",
                biom.Table(
                    inventory.values,
                    observation_ids=inventory.index,
                    sample_ids=inventory.columns,
                ),
            ),
            taxonomy=qiime2.Artifact.import_data(
                "FeatureData[Taxonomy]", self.taxonomy
            ),
            taxon_to_contig_map=qiime2.Artifact.import_data(
                "FeatureMap[TaxonomyToContigs]", self.mapping
            ),
        )
        self.assertEqual(
            str(result.feature_load.type), "FeatureTable[Frequency % Properties('tfa')]"
        )
        self.assertEqual(str(result.gene_taxonomy.type), "FeatureData[GeneTaxonomy]")
        table, mapping = self.estimate(abundance=self.abundance * 0.01)
        np.testing.assert_array_equal(
            result.feature_load.view(biom.Table).matrix_data.toarray(),
            table.matrix_data.toarray(),
        )
        pd.testing.assert_frame_equal(result.gene_taxonomy.view(pd.DataFrame), mapping)

    def test_invalid_abundances_rejected(self):
        for invalid in (-1, float("nan"), float("inf")):
            abundance = self.abundance.astype(float)
            abundance.loc["S1", "C1"] = invalid
            with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
                self.estimate(abundance=abundance)
