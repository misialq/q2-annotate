# -----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
"""Interactive exploration of sample-resolved TFA loads."""

import json
import shutil
from importlib import resources
from pathlib import Path

import biom
import numpy as np
import pandas as pd
from qiime2 import Metadata

from .types._format import _split_pair_id


def _visualization_data(feature_load: biom.Table, metadata: Metadata | None) -> dict:
    sample_ids = list(feature_load.ids(axis="sample"))
    matrix = feature_load.matrix_data.tocsr()
    pairs = []
    loads = []

    for row, observation_id in enumerate(feature_load.ids(axis="observation")):
        taxon_id, function_id = _split_pair_id(observation_id)
        start, end = matrix.indptr[row], matrix.indptr[row + 1]
        values = matrix.data[start:end]
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError("TFA loads must be finite and nonnegative.")
        pairs.append(
            {
                "taxon": taxon_id,
                "function": function_id,
                "total": float(values.sum()),
            }
        )
        loads.extend(
            [row, int(sample), float(value)]
            for sample, value in zip(matrix.indices[start:end], values)
            if value != 0
        )

    groups = {}
    if metadata is not None:
        frame = metadata.to_dataframe()
        for column, properties in metadata.columns.items():
            if properties.type == "categorical":
                series = frame[column].reindex(sample_ids)
                groups[column] = [
                    "(missing)" if pd.isna(value) else str(value) for value in series
                ]

    return {"samples": sample_ids, "pairs": pairs, "loads": loads, "groups": groups}


def explore_tfa(
    output_dir: str, feature_load: biom.Table, metadata: Metadata | None = None
) -> None:
    """Write an interactive TFA explorer with optional sample groups."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    assets = resources.files("q2_annotate") / "assets" / "tfa_explore"
    payload = json.dumps(
        _visualization_data(feature_load, metadata),
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    )
    # JSON inside a script element must not contain literal HTML delimiters.
    payload = (
        payload.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    )
    html = (assets / "index.html").read_text(encoding="utf-8")
    (output / "index.html").write_text(
        html.replace("__TFA_DATA__", payload), encoding="utf-8"
    )
    for filename in ("explore.js", "style.css"):
        shutil.copyfile(assets / filename, output / filename)
