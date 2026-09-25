# -----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import json
import re
from pathlib import Path

import biom
import pandas as pd
import qiime2
import scipy.sparse as sp
from qiime2.plugin.testing import TestPluginBase

from q2_annotate.plugin_setup import plugin
from q2_annotate.tfa.types._format import _pair_id
from q2_annotate.tfa.visualizer import explore_tfa


class TestTFAVisualizer(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def _table(self):
        data = pd.read_csv(
            self.get_data_path("tfa-table.tsv"),
            sep="\t",
            dtype={"taxon_id": str, "function_id": str},
        )
        samples = list(data.columns[2:])
        return biom.Table(
            sp.csr_matrix(data[samples].to_numpy(dtype=float)),
            observation_ids=[
                _pair_id(taxon, function)
                for taxon, function in zip(data.taxon_id, data.function_id)
            ],
            sample_ids=samples,
        )

    def _payload(self, html):
        match = re.search(
            r'<script id="tfa-data" type="application/json">(.*?)</script>',
            html,
            re.DOTALL,
        )
        self.assertIsNotNone(match)
        return json.loads(match.group(1))

    def test_sparse_payload_and_metadata_groups(self):
        metadata = qiime2.Metadata.load(self.get_data_path("sample-metadata.tsv"))
        explore_tfa(self.temp_dir.name, self._table(), metadata=metadata)
        output = Path(self.temp_dir.name)
        html = (output / "index.html").read_text(encoding="utf-8")
        payload = self._payload(html)
        self.assertEqual(payload["samples"], ["S1", "S2", "S3"])
        self.assertEqual([pair["total"] for pair in payload["pairs"]], [8, 4.25, 16])
        self.assertEqual(len(payload["loads"]), 6)
        self.assertEqual(payload["loads"][0], [0, 0, 5.0])
        self.assertEqual(
            payload["groups"],
            {
                "group": ["case", "control", "case"],
                "batch": ["B1", "B2", "B1"],
            },
        )
        self.assertIn("Explore taxon", html)
        self.assertTrue((output / "explore.js").is_file())
        self.assertTrue((output / "style.css").is_file())

    def test_optional_metadata_and_html_escaping(self):
        table = biom.Table(
            sp.csr_matrix([[2]]),
            observation_ids=[_pair_id("T<script>", "gene")],
            sample_ids=["S1"],
        )
        explore_tfa(self.temp_dir.name, table)
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        self.assertNotIn("T<script>", html)
        self.assertEqual(self._payload(html)["pairs"][0]["taxon"], "T<script>")
        self.assertEqual(self._payload(html)["groups"], {})

    def test_samples_missing_from_metadata_remain_grouped(self):
        metadata = qiime2.Metadata.load(
            self.get_data_path("sample-metadata-partial.tsv")
        )
        explore_tfa(self.temp_dir.name, self._table(), metadata=metadata)
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        self.assertEqual(
            self._payload(html)["groups"]["group"],
            ["case", "control", "(missing)"],
        )

    def test_action_registration_and_execution(self):
        self.assertIn("explore_tfa", plugin.visualizers)
        artifact = qiime2.Artifact.import_data("FeatureTable[TFA]", self._table())
        metadata = qiime2.Metadata.load(self.get_data_path("sample-metadata.tsv"))
        result = plugin.visualizers["explore_tfa"](
            feature_load=artifact, metadata=metadata
        )
        self.assertIsInstance(result.visualization, qiime2.Visualization)
        without_metadata = plugin.visualizers["explore_tfa"](feature_load=artifact)
        self.assertIsInstance(without_metadata.visualization, qiime2.Visualization)

    def test_taxonomy_labels_and_missing_assignments(self):
        taxonomy = pd.read_csv(
            self.get_data_path("taxonomy.tsv"), sep="\t", index_col=0
        )
        explore_tfa(self.temp_dir.name, self._table(), taxonomy=taxonomy)
        html = (Path(self.temp_dir.name) / "index.html").read_text(encoding="utf-8")
        pairs = self._payload(html)["pairs"]
        self.assertEqual(pairs[0]["taxon_id"], "T1")
        self.assertEqual(pairs[0]["taxon"], "d__Bacteria; g__Bacteroides; s__fragilis")
        self.assertEqual(pairs[0]["taxon_short"], "s__fragilis")
        self.assertEqual(pairs[2]["taxon"], "T2")
        self.assertEqual(pairs[2]["taxon_short"], "T2")

    def test_action_accepts_taxonomy_artifact(self):
        table = qiime2.Artifact.import_data("FeatureTable[TFA]", self._table())
        taxonomy = pd.read_csv(
            self.get_data_path("taxonomy.tsv"), sep="\t", index_col=0
        )
        taxonomy.index.name = "Feature ID"
        assignment = qiime2.Artifact.import_data("FeatureData[Taxonomy]", taxonomy)
        result = plugin.visualizers["explore_tfa"](
            feature_load=table, taxonomy=assignment
        )
        self.assertIsInstance(result.visualization, qiime2.Visualization)
