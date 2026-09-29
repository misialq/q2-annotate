# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from collections import defaultdict

import biom
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .types._format import GENE_TAXONOMY_COLUMNS
from .types._utils import _gene_taxonomy_id


def estimate_tfa(
    abundance_matrix: biom.Table,
    feature_inventory: biom.Table,
    taxonomy: pd.DataFrame,
    taxon_to_contig_map: dict,
) -> (biom.Table, pd.DataFrame):
    """Estimate sparse, sample-resolved taxon/function loads."""
    for table in (abundance_matrix, feature_inventory):
        data = table.matrix_data.data
        if not np.all(np.isfinite(data)) or np.any(data < 0):
            raise ValueError(
                "Abundances and gene counts must be finite and nonnegative."
            )
    # Build a mapping from contig_id to taxon_id using taxon_to_contig_map and taxonomy
    contig_to_taxon_id = {}
    for taxon_id, contigs in taxon_to_contig_map.items():
        if taxon_id in taxonomy.index:
            taxon_name = taxonomy.loc[taxon_id, "Taxon"]
            if pd.notnull(taxon_name):
                for cid in contigs:
                    if (
                        cid in contig_to_taxon_id
                        and contig_to_taxon_id[cid] != taxon_id
                    ):
                        raise ValueError(
                            f"Contig {cid!r} is assigned to multiple taxa."
                        )
                    contig_to_taxon_id[cid] = taxon_id

    # Find the intersection of contig IDs across all inputs
    common_contigs = (
        set(abundance_matrix.ids(axis="observation"))
        & set(feature_inventory.ids(axis="sample"))
        & set(contig_to_taxon_id.keys())
    )

    feature_ids = feature_inventory.ids(axis="observation")
    sample_ids = abundance_matrix.ids(axis="sample")

    # Sort common contigs to ensure deterministic ordering
    common_contigs = sorted(list(common_contigs))

    # Slice and align the abundance matrix (observations x samples)
    abund_contig_index = {
        cid: idx for idx, cid in enumerate(abundance_matrix.ids(axis="observation"))
    }
    abund_indices = [abund_contig_index[cid] for cid in common_contigs]
    abundance_matrix_sparse = abundance_matrix.matrix_data.tocsr()

    A = abundance_matrix_sparse[abund_indices, :]

    # Slice and align the feature inventory matrix. Its observations are
    # functions and its samples are contigs.
    feat_contig_index = {
        cid: idx for idx, cid in enumerate(feature_inventory.ids(axis="sample"))
    }
    feat_indices = [feat_contig_index[cid] for cid in common_contigs]
    feature_inventory_sparse = feature_inventory.matrix_data.tocsc()

    inventory = feature_inventory_sparse[:, feat_indices].T.tocsr()
    inventory.eliminate_zeros()

    # For each taxon, multiply its sparse function-by-contig inventory by
    # contig-by-sample abundances. Only observed taxon/function pairs become
    # rows; no dense taxon x function x sample cube is constructed.
    contig_rows_by_taxon = defaultdict(list)
    for row, contig_id in enumerate(common_contigs):
        contig_rows_by_taxon[contig_to_taxon_id[contig_id]].append(row)
    blocks = []
    pair_ids = []
    gene_taxonomy_rows = []
    for taxon_id, contig_rows in sorted(contig_rows_by_taxon.items()):
        taxon_inventory = inventory[contig_rows, :]
        feature_columns = np.unique(taxon_inventory.indices)
        if not len(feature_columns):
            continue
        loads = (taxon_inventory[:, feature_columns].T @ A[contig_rows, :]).tocsr()
        loads.eliminate_zeros()
        blocks.append(loads)
        for column in feature_columns:
            gene_id = feature_ids[column]
            pair_ids.append(_gene_taxonomy_id(taxon_id, gene_id))
            gene_taxonomy_rows.append(
                [taxon_id, gene_id, taxonomy.at[taxon_id, "Taxon"]]
            )

    result_matrix = (
        sp.vstack(blocks, format="csr")
        if blocks
        else sp.csr_matrix((0, len(sample_ids)))
    )
    if not np.all(np.isfinite(result_matrix.data)):
        raise ValueError("Estimated loads must be finite.")
    gene_taxonomy = pd.DataFrame(
        gene_taxonomy_rows,
        index=pd.Index(pair_ids, name="Feature ID"),
        columns=GENE_TAXONOMY_COLUMNS,
    )
    return (
        biom.Table(result_matrix, observation_ids=pair_ids, sample_ids=sample_ids),
        gene_taxonomy,
    )
