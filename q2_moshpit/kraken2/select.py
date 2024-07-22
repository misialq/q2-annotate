# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
from collections import deque
from typing import List

from q2_moshpit.kraken2.utils import (
    _find_lca, _taxon_to_list, _join_ranks
)
from q2_types.kraken2 import (
    Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat
)

import pandas as pd
import skbio

RANKS = 'dkpcofgs'


def _fill_unclassified(row):
    if row.isnull().all():
        row[0] = 'Unclassified'
    return row


def _find_lcas(taxa_list: List[pd.DataFrame], mode: str):
    """Find the least common ancestor in every DataFrame of taxa.

    Args:
        taxa_list (List[pd.DataFrame]): A list of taxonomy DataFrames.
        mode (str): The mode used to determine the least common ancestor.

    Returns:
        pd.DataFrame: A DataFrame containing the LCA of each feature (MAG).
    """
    methods = {
        'lca': _find_lca,
        # 'super': _find_super_lca,
        # 'majority': _find_lca_majority
    }
    func = methods[mode]
    taxa = pd.concat(taxa_list)

    # Convert taxonomies to list; optionally remove rank handle
    taxa['Taxon'] = taxa['Taxon'].apply(
        lambda x: _taxon_to_list(x, rank_handle=f'^[{RANKS[:-1]}]__|s1?__')
    )

    # Find LCA for every MAG
    results = {}
    for mag_id in taxa['mag_id'].unique():
        data = taxa[taxa['mag_id'] == mag_id]['Taxon']
        result = func(data)
        results[mag_id] = result

    results = pd.DataFrame.from_dict(results, orient='index')
    results = results.apply(_fill_unclassified, axis=1)
    results = results.apply(lambda x: x.tolist(), axis=1).to_frame()
    results.columns = ['Taxon']

    # Join ranks
    ranks = [*[f'{r}__' for r in RANKS], 'ssp__']
    results['Taxon'] = results['Taxon'].apply(
        lambda x: _join_ranks(x, ranks)
    )

    results.index.name = 'Feature ID'
    return results


def _add_unclassified_mags(
        table: pd.DataFrame,
        taxonomy: pd.DataFrame,
        reports: Kraken2ReportDirectoryFormat
) -> (pd.DataFrame, pd.Series):
    # identify which MAGs are entirely unclassified - they will be missing
    # from the table
    samples = [
        os.path.basename(f).replace(".report.txt", "")
        for f in glob.glob(os.path.join(reports.path, "*.report.txt"))
    ]
    unclassified = set(samples) - set(table.index)

    for mag in unclassified:
        report_fp = os.path.join(reports.path, f"{mag}.report.txt")
        with open(report_fp, "r") as f:
            line = f.readline()
            if "unclassified" in line:
                fraction = float(line.split("\t")[0])
                if fraction != 100.0:
                    raise ValueError(
                        f"Unclassified fraction for MAG '{mag}' is not "
                        f"100.0%: {fraction:.2f}%."
                    )
                taxonomy.loc[mag, "Taxon"] = "d__Unclassified"
            else:
                raise ValueError(
                    f"Unclassified line for MAG '{mag}' is missing from "
                    f"the Kraken 2 report."
                )
    return taxonomy


def kraken2_to_mag_features(
        reports: Kraken2ReportDirectoryFormat,
        hits: Kraken2OutputDirectoryFormat,
        coverage_threshold: float = 0.1,
        # lca_mode: str = 'lca'
) -> pd.DataFrame:
    table, taxonomy = kraken2_to_features(reports, coverage_threshold)

    taxa_list = []
    # convert IDs to match MAGs instead of taxids/db ids
    for mag_id in table.index:
        kraken_table_fp = (hits.path / f'{mag_id}.output.txt')
        hits_df = pd.read_csv(
            kraken_table_fp, sep='\t', header=None, dtype='str'
        )
        MAG_COL = 1
        TAXA_COL = 2

        mag_series = table.loc[mag_id, :]
        mag_obs = mag_series[mag_series != 0]
        merged_df = hits_df.join(mag_obs, on=TAXA_COL, how='right')
        merged_df = merged_df.join(taxonomy, on=TAXA_COL, how='left')

        new_taxa = merged_df[[MAG_COL, 'Taxon']].set_index(MAG_COL)
        new_taxa.index.name = 'Feature ID'
        new_taxa['mag_id'] = mag_id
        taxa_list.append(new_taxa)

    taxonomy = _find_lcas(taxa_list, mode='lca')
    return _add_unclassified_mags(table, taxonomy, reports)


def kraken2_to_features(reports: Kraken2ReportDirectoryFormat,
                        coverage_threshold: float = 0.1) \
        -> (pd.DataFrame, pd.DataFrame):

    rows = []
    trees = []
    unclassified = {}
    for relpath, df in reports.reports.iter_views(pd.DataFrame):
        sample_id = os.path.basename(relpath).replace(".report.txt", "")

        filtered = df[df['perc_frags_covered'] >= coverage_threshold]
        unclassified[sample_id] = filtered[filtered['name'] == 'unclassified']
        tree = _kraken_to_ncbi_tree(filtered)
        tips = _ncbi_tree_to_tips(tree)
        if tips:
            table_row = pd.Series(True, index=_ncbi_tree_to_tips(tree))
            table_row.name = sample_id
            rows.append(table_row)
        trees.append(tree)

    full_tree = _combine_ncbi_trees(trees)

    table = pd.DataFrame(rows).fillna(False)
    taxonomy = _to_taxonomy(full_tree)
    # filter taxonomy to only IDs in table
    # use list to avoid index name change
    taxonomy = taxonomy.loc[list(table.columns)]

    return table, taxonomy


def _get_indentation(string, indent=2):
    return (len(string) - len(string.lstrip(' '))) // indent


def _kraken_to_ncbi_tree(df):
    tree = skbio.TreeNode()
    stack = deque([(0, tree)])
    for _, row in df.iterrows():
        r = row['rank']
        label = row['name']
        otu = str(row['taxon_id'])

        if r in ('U', 'R'):
            continue  # unclassified or root

        indent = _get_indentation(label)
        name = f"{r.lower()}__{label.strip()}"
        node = skbio.TreeNode(name=name, length=0.0)

        # Don't include internal non-strain infra-clades as tips
        if len(r) == 1 or r.startswith('S'):
            id_node = skbio.TreeNode(name=otu, length=0.0)
            node.length = 1.0  # not infra-clade, so give it a length
            node.append(id_node)

        parent_indent, parent_node = stack[-1]
        if parent_indent >= indent and parent_node.children:
            parent_node.children[0].is_actual_tip = True

        while parent_indent >= indent:
            stack.pop()
            parent_indent, parent_node = stack[-1]

        parent_node.append(node)
        stack.append((indent, node))

    # last entry is always a tip
    _, parent_node = stack[-1]
    # It is possible for the last row to be an infra-clade tip,
    # so walk backwards up the stack until a standard node (with length)
    # is found
    while stack and parent_node.length == 0:
        _, parent_node = stack.pop()

    if parent_node.children:
        parent_node.children[0].is_actual_tip = True

    return tree


def _combine_ncbi_trees(trees):
    full_tree = trees[0]
    for tree in trees[1:]:
        for tip in list(tree.tips()):
            try:
                # check if taxid is already in this tree
                full_tree.find(tip.name)
                continue  # for clarity
            except skbio.tree.MissingNodeError:
                parents = list(tip.ancestors())[:-1]  # ignore unnamed root
                matching = full_tree
                while parents:
                    node = parents.pop()
                    try:
                        matching = matching.find(node.name)
                    except skbio.tree.MissingNodeError:
                        matching.append(node)
                        break
    return full_tree


def _ncbi_tree_to_tips(tree):
    return [n.name for n in tree.tips() if hasattr(n, 'is_actual_tip')]


def _pad_ranks(ranks):
    order = ['d', 'k', 'p', 'c', 'o', 'f', 'g', 's']
    available = {}
    taxonomy = []

    for rank in reversed(ranks):
        r, label = rank.split('__', 1)
        available[r] = label
        if len(r) > 1 and r.startswith('s'):
            taxonomy.append(rank)

    for r in reversed(order):
        label = available.get(r)

        if label is not None:
            taxonomy.append(f'{r}__{label}')
            last_good_label = f'{r}__{label}'
        elif taxonomy:
            if r == 'k' and available.get('d') in ('Bacteria', 'Archaea'):
                # it is too strange to propagate upwards for these 'kingdoms'
                taxonomy.append(f"k__{available['d']}")
            else:
                # smear the most specific label we have upwards
                taxonomy.append(f'{r}__containing {last_good_label}')

    return ';'.join(reversed(taxonomy))


def _to_taxonomy(tree):
    rows = [(node.name, _pad_ranks(ranks))
            for node, ranks in tree.to_taxonomy()]
    taxonomy = pd.DataFrame(rows, columns=['Feature ID', 'Taxon'])
    taxonomy = taxonomy.set_index('Feature ID')

    return taxonomy
