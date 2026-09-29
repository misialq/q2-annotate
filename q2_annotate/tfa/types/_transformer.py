# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd

from q2_annotate.plugin_setup import plugin

from ._format import GeneTaxonomyFormat
from ._utils import _validate_gene_taxonomy


@plugin.register_transformer
def _gene_taxonomy_to_format(frame: pd.DataFrame) -> GeneTaxonomyFormat:
    _validate_gene_taxonomy(frame)
    result = GeneTaxonomyFormat()
    with result.open() as handle:
        frame.to_csv(handle, sep="\t", index=True, index_label="Feature ID")
    return result


@plugin.register_transformer
def _gene_taxonomy_to_dataframe(source: GeneTaxonomyFormat) -> pd.DataFrame:
    with source.open() as handle:
        return pd.read_csv(
            handle, sep="\t", index_col=0, dtype=str, keep_default_na=False
        )
