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

    def test_tfa(self):
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
