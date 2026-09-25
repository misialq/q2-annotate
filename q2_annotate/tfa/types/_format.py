# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json

import h5py
import numpy as np
from qiime2.core.exceptions import ValidationError
from qiime2.plugin import model
from q2_types.feature_table import BIOMV210Format


def _pair_id(taxon_id: str, feature_id: str) -> str:
    """Encode both feature dimensions without assuming anything about IDs."""
    return json.dumps(
        [str(taxon_id), str(feature_id)], ensure_ascii=False, separators=(",", ":")
    )


def _split_pair_id(observation_id: str) -> tuple[str, str]:
    try:
        pair = json.loads(observation_id)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid TFA feature ID: {observation_id!r}") from error
    if (
        not isinstance(pair, list)
        or len(pair) != 2
        or not all(isinstance(value, str) for value in pair)
    ):
        raise ValueError(f"Invalid TFA feature ID: {observation_id!r}")
    return pair[0], pair[1]


class TFAFeatureTableFormat(model.BinaryFileFormat):
    """BIOM v2.1 with encoded taxon/function feature IDs."""

    def open(self):
        return h5py.File(str(self), mode=self._mode)

    def _validate_(self, level):
        try:
            with h5py.File(str(self), mode="r") as handle:
                for group in BIOMV210Format.groups:
                    if group not in handle:
                        raise ValidationError(f"Missing BIOM group: {group}")
                for dataset in BIOMV210Format.datasets:
                    if dataset not in handle:
                        raise ValidationError(f"Missing BIOM dataset: {dataset}")
                for attribute in BIOMV210Format.attrs:
                    if attribute not in handle.attrs:
                        raise ValidationError(f"Missing BIOM attribute: {attribute}")
                ids = handle["observation/ids"]
                values = handle["observation/matrix/data"]
                id_count = min(len(ids), 100) if level == "min" else len(ids)
                value_count = min(len(values), 100) if level == "min" else len(values)
                for start in range(0, id_count, 8192):
                    for raw_id in ids[start : min(start + 8192, id_count)]:
                        try:
                            _split_pair_id(raw_id.decode("utf-8"))
                        except (UnicodeDecodeError, ValueError) as error:
                            raise ValidationError(
                                f"Invalid taxon/function pair ID: {raw_id!r}"
                            ) from error
                for start in range(0, value_count, 8192):
                    data = values[start : min(start + 8192, value_count)]
                    if not np.all(np.isfinite(data)) or np.any(data < 0):
                        raise ValidationError(
                            "TFA loads must be finite and nonnegative."
                        )
        except OSError as error:
            raise ValidationError("Expected a BIOM v2.1 TFA feature table.") from error


TFAFeatureTableDirFmt = model.SingleFileDirectoryFormat(
    "TFAFeatureTableDirFmt",
    "tfa-table.biom",
    TFAFeatureTableFormat,
)
