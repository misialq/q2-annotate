# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from typing import Union

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt

from q2_annotate._utils import run_command
from .types import GUNCResultsDirectoryFormat, GUNCDatabaseDirFmt


def download_gunc_db(database: str = "progenomes") -> GUNCDatabaseDirFmt:
    """Download the GUNC reference database."""

    db = GUNCDatabaseDirFmt()
    cmd = ["gunc", "download_db", str(db.path)]
    if database:
        cmd.extend(["--database", database])

    run_command(cmd, verbose=True)
    return db


def run_gunc(
    mags: Union[MAGSequencesDirFmt, MultiMAGSequencesDirFmt],
    db: GUNCDatabaseDirFmt,
    threads: int = 1,
    sensitive: bool = False,
    detailed_output: bool = False,
    contig_taxonomy_output: bool = False,
    use_species_level: bool = False,
    min_mapped_genes: int = 11,
) -> GUNCResultsDirectoryFormat:
    """Run GUNC on the provided MAGs."""

    results = GUNCResultsDirectoryFormat()

    cmd = [
        "gunc",
        "run",
        "--input_dir",
        str(mags.path),
        "--db_file",
        str(db.path),
        "--out_dir",
        str(results.path),
        "--threads",
        str(threads),
    ]

    if sensitive:
        cmd.append("--sensitive")
    if detailed_output:
        cmd.append("--detailed_output")
    if contig_taxonomy_output:
        cmd.append("--contig_taxonomy_output")
    if use_species_level:
        cmd.append("--use_species_level")
    if min_mapped_genes is not None:
        cmd.extend(["--min_mapped_genes", str(min_mapped_genes)])

    run_command(cmd, verbose=True)
    return results
