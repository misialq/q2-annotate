# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import hashlib
import json

import pandas as pd

GENE_TAXONOMY_COLUMNS = ["Taxon ID", "Gene ID", "Taxon"]


def _gene_taxonomy_id(taxon_id: str, gene_id: str) -> str:
    """Identify a pair independently of its order or surrounding dataset."""
    pair = json.dumps(
        [str(taxon_id), str(gene_id)], ensure_ascii=False, separators=(",", ":")
    )
    return "gt_" + hashlib.sha256(pair.encode("utf-8")).hexdigest()


def _validate_gene_taxonomy(frame: pd.DataFrame) -> None:
    if list(frame.columns) != GENE_TAXONOMY_COLUMNS:
        raise ValueError("Gene taxonomy columns must be Taxon ID, Gene ID, Taxon.")
    if not frame.index.is_unique:
        raise ValueError("Gene taxonomy feature IDs must be unique.")
    for values in (frame.index, frame["Taxon ID"], frame["Gene ID"]):
        if any(not isinstance(value, str) or not value.strip() for value in values):
            raise ValueError("Gene taxonomy IDs must be nonempty strings.")
    if frame.duplicated(subset=["Taxon ID", "Gene ID"]).any():
        raise ValueError("Gene taxonomy taxon/gene pairs must be unique.")
    if any(not isinstance(value, str) for value in frame["Taxon"]):
        raise ValueError("Gene taxonomy labels must be strings (empty is allowed).")
