# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import biom
import pandas as pd
from rachis.plugin.testing import TestPluginBase

from q2_annotate.tfa import estimate_tfa


class TestAbundance(TestPluginBase):
    package = "q2_annotate.abundance.tests"

    def setUp(self):
        super().setUp()
        # df_tpm: Sample x Contig (Abundance)
        self.df_tpm = pd.DataFrame(
            {
                "C1": [100, 0, 0],
                "C2": [50, 0, 0],
                "C3": [0, 120, 0],
                "C4": [0, 40, 0],
                "C5": [0, 0, 90],
                "C6": [0, 0, 60],
            },
            index=pd.Index(["S1", "S2", "S3"], name="id"),
        )

        # df_amr: Contig x AMR (Gene counts)
        self.df_amr = pd.DataFrame(
            {
                "bla_TEM": [2, 1, 2, 1, 2, 1],
                "vanA": [1, 0, 1, 0, 1, 0],
                "mecA": [0, 3, 0, 3, 0, 3],
            },
            index=pd.Index(["C1", "C2", "C3", "C4", "C5", "C6"], name="id"),
        )

        # df_taxon: Contig x Taxon (Mapping)
        self.df_taxon = pd.DataFrame(
            {
                "Taxon": [
                    "Escherichia coli",
                    "Klebsiella pneumoniae",
                ]
            },
            index=pd.Index(["T1", "T2"], name="id"),
        )

        self.taxon_to_contig_map = {
            "T1": ["C1", "C3", "C5"],
            "T2": ["C2", "C4", "C6"],
        }

    def _estimate_from_frames(self, abundance=None, taxonomy=None, mapping=None):
        if abundance is None:
            abundance = self.df_tpm
        if taxonomy is None:
            taxonomy = self.df_taxon
        if mapping is None:
            mapping = self.taxon_to_contig_map

        abundance = abundance.T
        abundance_biom = biom.Table(
            abundance.values,
            observation_ids=list(abundance.index),
            sample_ids=list(abundance.columns),
        )
        inventory = self.df_amr.T
        inventory_biom = biom.Table(
            inventory.values,
            observation_ids=list(inventory.index),
            sample_ids=list(inventory.columns),
        )
        return estimate_tfa(abundance_biom, inventory_biom, taxonomy, mapping)

    def test_tfa(self):
        """Aggregate abundance-weighted feature counts by taxon."""
        # Convert df_tpm to biom.Table
        # (transposed so that Contigs are observations, Samples are samples)
        tpm_df = self.df_tpm.T
        tpm_biom = biom.Table(
            tpm_df.values,
            observation_ids=list(tpm_df.index),
            sample_ids=list(tpm_df.columns),
        )

        # Convert df_amr to biom.Table (AMR genes as observations, Contigs as samples)
        amr_df = self.df_amr.T
        amr_biom = biom.Table(
            amr_df.values,
            observation_ids=list(amr_df.index),
            sample_ids=list(amr_df.columns),
        )

        # Construct the corresponding taxonomy DataFrame for self.taxon_to_contig_map
        # T1 represents "Escherichia coli", T2 represents "Klebsiella pneumoniae"
        taxonomy_df = pd.DataFrame(
            {"Taxon": ["Escherichia coli", "Klebsiella pneumoniae"]},
            index=pd.Index(["T1", "T2"], name="id"),
        )

        obs_biom = estimate_tfa(
            tpm_biom, amr_biom, taxonomy_df, self.taxon_to_contig_map
        )
        obs_df = obs_biom.to_dataframe(dense=True)
        obs_df.index = obs_df.index.map(taxonomy_df["Taxon"])

        exp = pd.DataFrame(
            {"bla_TEM": [620, 150], "mecA": [0, 450], "vanA": [310, 0]},
            index=pd.Index(["Escherichia coli", "Klebsiella pneumoniae"], name="Taxon"),
        )
        exp.columns.name = "feature_id"

        obs_df.index.name = "Taxon"
        obs_df.columns.name = "feature_id"
        obs_df = obs_df.reindex(index=exp.index, columns=exp.columns)

        pd.testing.assert_frame_equal(obs_df, exp, check_dtype=False)

    def test_tfa_no_common_contigs(self):
        """Return no taxa but retain feature IDs when no contigs overlap."""
        tpm_df = self.df_tpm.T
        tpm_biom = biom.Table(
            tpm_df.values,
            observation_ids=list(tpm_df.index),
            sample_ids=list(tpm_df.columns),
        )

        amr_df = self.df_amr.T
        amr_biom = biom.Table(
            amr_df.values,
            observation_ids=list(amr_df.index),
            sample_ids=list(amr_df.columns),
        )

        obs_biom = estimate_tfa(
            tpm_biom,
            amr_biom,
            self.df_taxon,
            {"T1": ["not-a-contig"]},
        )

        self.assertEqual(obs_biom.shape, (0, len(amr_df.index)))
        self.assertEqual(list(obs_biom.ids(axis="observation")), [])
        self.assertEqual(list(obs_biom.ids(axis="sample")), list(amr_df.index))

    def test_tfa_ignores_extra_contigs_in_abundance_and_taxonomy(self):
        """Ignore contigs outside the three-way input intersection."""
        abundance = self.df_tpm.assign(C7=[1000, 1000, 1000])
        taxonomy = pd.DataFrame(
            {
                "Taxon": [
                    "Escherichia coli",
                    "Klebsiella pneumoniae",
                    "Extra taxon",
                ]
            },
            index=["T1", "T2", "T3"],
        )
        mapping = {
            "T1": ["C1", "C3", "C5"],
            "T2": ["C2", "C4", "C6"],
            "T3": ["C7", "C8"],
        }

        observed = self._estimate_from_frames(
            abundance=abundance, taxonomy=taxonomy, mapping=mapping
        ).to_dataframe(dense=True)
        expected = pd.DataFrame(
            {
                "bla_TEM": [620, 150],
                "vanA": [310, 0],
                "mecA": [0, 450],
            },
            index=["T1", "T2"],
        )

        pd.testing.assert_frame_equal(observed, expected, check_dtype=False)

    def test_tfa_ignores_abundance_sample_ids(self):
        """Compute the same result when abundance sample labels differ."""
        abundance = self.df_tpm.copy()
        abundance.index = ["other-1", "other-2", "other-3"]

        observed = self._estimate_from_frames(
            abundance=abundance
        ).to_dataframe(dense=True)
        expected = pd.DataFrame(
            {
                "bla_TEM": [620, 150],
                "vanA": [310, 0],
                "mecA": [0, 450],
            },
            index=["T1", "T2"],
        )

        pd.testing.assert_frame_equal(observed, expected, check_dtype=False)

    def test_tfa_keeps_zero_abundance_contig(self):
        """Retain a taxon whose only contig has zero abundance."""
        abundance = pd.DataFrame(
            {"C1": [0, 0, 0], "C2": [50, 0, 0]},
            index=["S1", "S2", "S3"],
        )

        observed = self._estimate_from_frames(
            abundance=abundance, mapping={"T1": ["C1"], "T2": ["C2"]}
        ).to_dataframe(dense=True)
        expected = pd.DataFrame(
            {
                "bla_TEM": [0, 50],
                "vanA": [0, 0],
                "mecA": [0, 150],
            },
            index=["T1", "T2"],
        )

        pd.testing.assert_frame_equal(observed, expected, check_dtype=False)

    def test_tfa_all_zero_abundance(self):
        """Retain mapped taxa and features when all abundances are zero."""
        abundance = pd.DataFrame(
            {
                "C1": [0, 0, 0],
                "C2": [0, 0, 0],
                "C3": [0, 0, 0],
                "C4": [0, 0, 0],
                "C5": [0, 0, 0],
                "C6": [0, 0, 0],
            },
            index=["S1", "S2", "S3"],
        )

        observed = self._estimate_from_frames(abundance=abundance)
        expected = pd.DataFrame(
            {
                "bla_TEM": [0, 0],
                "vanA": [0, 0],
                "mecA": [0, 0],
            },
            index=["T1", "T2"],
        )

        pd.testing.assert_frame_equal(
            observed.to_dataframe(dense=True), expected, check_dtype=False
        )
        self.assertEqual(observed.matrix_data.nnz, 0)

    def test_tfa_preserves_extreme_abundances(self):
        """Preserve tiny and large abundance contributions."""
        abundance = pd.DataFrame(
            {
                "C1": [1e-9, 0, 0],
                "C2": [1e12, 0, 0],
                "C3": [0, 0, 0],
                "C4": [0, 0, 0],
                "C5": [0, 0, 0],
                "C6": [0, 0, 0],
            },
            index=["S1", "S2", "S3"],
        )

        observed = self._estimate_from_frames(
            abundance=abundance
        ).to_dataframe(dense=True)
        expected = pd.DataFrame(
            {
                "bla_TEM": [2e-9, 1e12],
                "vanA": [1e-9, 0],
                "mecA": [0, 3e12],
            },
            index=["T1", "T2"],
        )

        pd.testing.assert_frame_equal(
            observed, expected, check_exact=False, rtol=1e-12, atol=0
        )

    def test_tfa_singleton_taxa(self):
        """Keep each contig's feature load in its own taxon."""
        taxonomy = pd.DataFrame(
            {
                "Taxon": [
                    "Taxon 1",
                    "Taxon 2",
                    "Taxon 3",
                    "Taxon 4",
                    "Taxon 5",
                    "Taxon 6",
                ]
            },
            index=["T1", "T2", "T3", "T4", "T5", "T6"],
        )
        mapping = {
            "T1": ["C1"],
            "T2": ["C2"],
            "T3": ["C3"],
            "T4": ["C4"],
            "T5": ["C5"],
            "T6": ["C6"],
        }

        observed = self._estimate_from_frames(
            taxonomy=taxonomy,
            mapping=mapping,
        ).to_dataframe(dense=True)
        expected = pd.DataFrame(
            {
                "bla_TEM": [200, 50, 240, 40, 180, 60],
                "vanA": [100, 0, 120, 0, 90, 0],
                "mecA": [0, 150, 0, 120, 0, 180],
            },
            index=["T1", "T2", "T3", "T4", "T5", "T6"],
        )

        pd.testing.assert_frame_equal(observed, expected, check_dtype=False)

    def test_tfa_super_taxon(self):
        """Sum every contig's feature load into one taxon."""
        observed = self._estimate_from_frames(
            mapping={"T1": ["C1", "C2", "C3", "C4", "C5", "C6"]}
        ).to_dataframe(dense=True)
        expected = pd.DataFrame(
            {"bla_TEM": [770], "vanA": [310], "mecA": [450]},
            index=["T1"],
        )

        pd.testing.assert_frame_equal(observed, expected, check_dtype=False)
