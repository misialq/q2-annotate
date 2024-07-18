# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess

import pandas as pd

from q2_types.feature_data import BLAST6
from q2_types.feature_data_mag import OrthologAnnotationDirFmt
from q2_types.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types.reference_db import EggnogRefDirFmt
from q2_types.sample_data import SampleData


def _annotate_seed_orthologs_runner(
        seed_ortholog, eggnog_db, sample_label, output_loc,
        db_in_memory, num_cpus
):
    # at this point instead of being able to specify the type of target
    # orthologs, we want to annotate _all_.

    cmds = ['emapper.py', '-m', 'no_search', '--annotate_hits_table',
            str(seed_ortholog), '--data_dir', str(eggnog_db),
            '-o', str(sample_label), '--output_dir', str(output_loc),
            '--cpu', str(num_cpus)]
    if db_in_memory:
        cmds.append('--dbmem')

    subprocess.run(cmds, check=True, cwd=str(output_loc))


def _eggnog_annotate(
        eggnog_hits: SeedOrthologDirFmt,
        eggnog_db: EggnogRefDirFmt,
        db_in_memory: bool = False,
        num_cpus: int = 1
) -> OrthologAnnotationDirFmt:

    eggnog_db_fp = eggnog_db.path

    result = OrthologAnnotationDirFmt()

    # run analysis
    for relpath, obj_path in eggnog_hits.seed_orthologs.iter_views(
            OrthologFileFmt):
        sample_label = str(relpath).rsplit(r'.', 2)[0]

        _annotate_seed_orthologs_runner(seed_ortholog=obj_path,
                                        eggnog_db=eggnog_db_fp,
                                        sample_label=sample_label,
                                        output_loc=result,
                                        db_in_memory=db_in_memory,
                                        num_cpus=num_cpus)

    return result


def eggnog_annotate(
        ctx,
        eggnog_hits,
        eggnog_db,
        db_in_memory=False,
        num_cpus=1,
        num_partitions=None
):
    _eggnog_annotate = ctx.get_action("moshpit", "_eggnog_annotate")
    collate_annotations = ctx.get_action("moshpit", "collate_annotations")

    if eggnog_hits.type <= SampleData[BLAST6]:
        partition_method = ctx.get_action("moshpit", "partition_orthologs")
    else:
        raise NotImplementedError()

    (partitioned_orthologs, ) = partition_method(eggnog_hits, num_partitions)

    annotations = []
    for orthologs in partitioned_orthologs.values():
        (annotation, ) = _eggnog_annotate(
            orthologs, eggnog_db, db_in_memory, num_cpus
        )
        annotations.append(annotation)

    (collated_annotations, ) = collate_annotations(annotations)
    return collated_annotations


def _extract_generic(
        data: pd.DataFrame, column: str, transform_func: callable
) -> pd.Series:
    data = data.set_index("seed_ortholog", inplace=False)[column]
    data = data.apply(transform_func if data.notna().any() else lambda x: x)
    data = data.stack().reset_index(level=1, drop=True)
    data = data.value_counts()
    for char in ("", "-"):
        if char in data.index:
            data = data.drop(char, axis=0, inplace=False)
    return data


def _extract_cog(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "COG_category", lambda x: pd.Series(list(x))
    )


def _extract_kegg_ko(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "KEGG_ko", lambda x: pd.Series([i[3:] for i in x.split(",")])
    )


def _extract_kegg_pathway(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "KEGG_Pathway", lambda x: pd.Series(
            [i for i in x.split(",") if i.startswith("map")]
        )
    )


def _extract_kegg_module(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "KEGG_Module", lambda x: pd.Series(x.split(","))
    )


def _extract_kegg_reaction(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "KEGG_Reaction", lambda x: pd.Series(x.split(","))
    )


def _extract_brite(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "BRITE", lambda x: pd.Series(x.split(","))
    )


def _extract_caz(data: pd.DataFrame) -> pd.Series:
    return _extract_generic(
        data, "CAZy", lambda x: pd.Series(x.split(","))
    )


def _filter(
        data: pd.DataFrame, max_evalue: float, min_score: float
) -> pd.DataFrame:
    data = data[(data["evalue"] <= max_evalue) & (data["score"] >= min_score)]
    if len(data) == 0:
        raise ValueError(
            "E-value/score filtering resulted in an empty table - "
            "please adjust your thresholds and try again."
        )
    return data


def extract_annotations(
        ortholog_annotations: OrthologAnnotationDirFmt,
        annotation: str,
        max_evalue: float = 1.0,
        min_score: float = 0.0
) -> pd.DataFrame:
    extract_method = globals().get(f"_extract_{annotation}")
    if not extract_method:
        raise NotImplementedError(
            f"Annotation method {annotation} not supported"
        )

    annotations = []
    for _id, fp in ortholog_annotations.annotation_dict().items():
        annot_df = pd.read_csv(
            fp, sep="\t", skiprows=4, index_col=0
        )  # skip the first 4 rows as they contain comments
        annot_df = annot_df.iloc[:-3, :]  # remove the last 3 comment rows
        annot_df = _filter(annot_df, max_evalue, min_score)
        annot_df = extract_method(annot_df)
        annot_df.name = _id
        annotations.append(annot_df)

    result = pd.concat(annotations, axis=1).fillna(0).T
    result.index.name = "id"
    return result
