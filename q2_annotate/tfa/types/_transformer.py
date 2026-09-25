# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import biom
import qiime2

from q2_annotate.plugin_setup import plugin

from ._format import TFAFeatureTableFormat


@plugin.register_transformer
def _biom_to_tfa_format(table: biom.Table) -> TFAFeatureTableFormat:
    result = TFAFeatureTableFormat()
    with result.open() as handle:
        table.to_hdf5(handle, generated_by=f"qiime2 {qiime2.__version__}")
    return result


@plugin.register_transformer
def _tfa_format_to_biom(source: TFAFeatureTableFormat) -> biom.Table:
    with source.open() as handle:
        table = biom.Table.from_hdf5(handle)
    # BIOM axis metadata is not part of this semantic type.
    return biom.Table(
        table.matrix_data,
        observation_ids=table.ids(axis="observation"),
        sample_ids=table.ids(axis="sample"),
    )
