# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from qiime2.plugin import Properties, ValidationError
from q2_types.feature_data import FeatureData, Taxonomy

from q2_annotate.plugin_setup import plugin
from . import TFA_TAXONOMY_COLUMNS


@plugin.register_validator(FeatureData[Taxonomy % Properties("tfa")])
def validate_tfa_taxonomy(data: pd.DataFrame, level):
    """Require a unique, labeled taxon/gene mapping for TFA feature IDs."""
    if not set(TFA_TAXONOMY_COLUMNS).issubset(data.columns):
        raise ValidationError(
            "TFA taxonomy requires Taxon, Taxon ID, and Gene ID columns."
        )

    # Minimal validation checks a prefix; full validation checks every pair.
    if level == "min":
        data = data.iloc[:100]

    for identifiers in (data.index.to_series(), data["Taxon ID"], data["Gene ID"]):
        invalid = identifiers.isna() | identifiers.astype(str).str.strip().eq("")
        if invalid.any():
            raise ValidationError(
                "TFA taxonomy feature IDs, taxon IDs, and gene IDs must be nonempty. "
                f"Invalid feature IDs: {data.index[invalid].tolist()}"
            )
    invalid = data["Taxon"].isna() | data["Taxon"].str.strip().eq("")
    if invalid.any():
        raise ValidationError(
            "TFA taxonomy taxon labels must be nonempty. "
            f"Invalid feature IDs: {data.index[invalid].tolist()}"
        )

    invalid = data.index.duplicated(keep=False) | data[
        ["Taxon ID", "Gene ID"]
    ].duplicated(keep=False)
    if invalid.any():
        raise ValidationError(
            "TFA taxonomy feature IDs and taxon/gene pairs must be unique."
            f" Duplicate feature IDs: {data.index[invalid].tolist()}"
        )
