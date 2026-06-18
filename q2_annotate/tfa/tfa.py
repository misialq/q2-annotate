# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import biom
import numpy as np
import pandas as pd
import scipy.sparse as sp


def estimate_tfa(
    abundance_matrix: biom.Table,
    feature_inventory: biom.Table,
    taxonomy: pd.DataFrame,
    taxon_to_contig_map: dict,
) -> biom.Table:
    # Build a mapping from contig_id to taxon_id using taxon_to_contig_map and taxonomy
    contig_to_taxon_id = {}
    for taxon_id, contigs in taxon_to_contig_map.items():
        if taxon_id in taxonomy.index:
            taxon_name = taxonomy.loc[taxon_id, "Taxon"]
            if pd.notnull(taxon_name):
                for cid in contigs:
                    contig_to_taxon_id[cid] = taxon_id

    # Find the intersection of contig IDs across all inputs
    common_contigs = set(abundance_matrix.ids(axis="observation")) & \
                     set(feature_inventory.ids(axis="sample")) & \
                     set(contig_to_taxon_id.keys())

    feature_ids = feature_inventory.ids(axis="observation")

    if not common_contigs:
        return biom.Table(
            sp.csr_matrix((0, len(feature_ids))),
            observation_ids=[],
            sample_ids=feature_ids,
        )

    # Sort common contigs to ensure deterministic ordering
    common_contigs = sorted(list(common_contigs))

    # Slice and align the abundance matrix (observations x samples)
    abund_contig_index = {
        cid: idx for idx, cid in enumerate(abundance_matrix.ids(axis="observation"))
    }
    abund_indices = [abund_contig_index[cid] for cid in common_contigs]

    if sp.issparse(abundance_matrix.matrix_data):
        abundance_matrix_sparse = abundance_matrix.matrix_data.tocsr()
    else:
        abundance_matrix_sparse = abundance_matrix.matrix_data

    A = abundance_matrix_sparse[abund_indices, :]

    # Slice and align the feature inventory matrix (observations x features, with contigs as samples)
    feat_contig_index = {
        cid: idx for idx, cid in enumerate(feature_inventory.ids(axis="sample"))
    }
    feat_indices = [feat_contig_index[cid] for cid in common_contigs]

    if sp.issparse(feature_inventory.matrix_data):
        feature_inventory_sparse = feature_inventory.matrix_data.tocsc()
    else:
        feature_inventory_sparse = feature_inventory.matrix_data

    I_sliced = feature_inventory_sparse[:, feat_indices]
    I = I_sliced.T

    # Compute total abundance of each contig across all samples (row sum of A)
    if sp.issparse(A):
        a_c = np.asarray(A.sum(axis=1)).flatten()
    else:
        a_c = A.sum(axis=1).flatten()

    # Scale the feature inventory by contig abundances
    if sp.issparse(I):
        M = sp.diags(a_c) @ I
    else:
        M = a_c[:, np.newaxis] * I

    # Group by taxon ID using a sparse grouping matrix G
    taxon_ids_list = np.array([contig_to_taxon_id[cid] for cid in common_contigs])
    unique_taxon_ids, taxon_indices = np.unique(taxon_ids_list, return_inverse=True)

    row_indices = taxon_indices
    col_indices = np.arange(len(common_contigs))
    data = np.ones(len(common_contigs))
    G = sp.coo_matrix(
        (data, (row_indices, col_indices)),
        shape=(len(unique_taxon_ids), len(common_contigs)),
    ).tocsr()

    # Compute total summed load per taxon ID
    result_matrix = G @ M

    # Return the resulting biom Table
    return biom.Table(
        result_matrix,
        observation_ids=list(unique_taxon_ids),
        sample_ids=list(feature_ids),
    )
