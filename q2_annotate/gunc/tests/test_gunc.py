# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase
from q2_annotate.gunc.gunc import run_gunc, download_gunc_db
from q2_annotate.gunc.types import GUNCResultsDirectoryFormat, GUNCDatabaseDirFmt
from q2_types.feature_data_mag import MAGSequencesDirFmt


class TestGUNC(TestPluginBase):
    package = "q2_annotate.gunc.tests"

    @patch("subprocess.run")
    def test_download_gunc_db(self, subp_run):
        db = download_gunc_db()
        subp_run.assert_called_once_with(
            ["gunc", "download_db", str(db.path), "--database", "progenomes"],
            check=True,
        )
        self.assertIsInstance(db, GUNCDatabaseDirFmt)

    @patch("subprocess.run")
    def test_run_gunc(self, subp_run):
        mags = MAGSequencesDirFmt(
            path=self.get_data_path("mags/dir_with_1_mag"), mode="r"
        )
        db = GUNCDatabaseDirFmt(path=self.get_data_path("db"), mode="r")
        results = run_gunc(
            mags=mags,
            db=db,
            threads=4,
            sensitive=True,
            detailed_output=True,
            contig_taxonomy_output=True,
            use_species_level=True,
            min_mapped_genes=42,
        )

        subp_run.assert_called_once_with(
            [
                "gunc",
                "run",
                "--input_dir",
                mags.path,
                "--db_file",
                str(db.path),
                "--out_dir",
                str(results.path),
                "--threads",
                "4",
                "--sensitive",
                "--detailed_output",
                "--contig_taxonomy_output",
                "--use_species_level",
                "--min_mapped_genes",
                "42",
            ],
            check=True,
        )
        self.assertIsInstance(results, GUNCResultsDirectoryFormat)
