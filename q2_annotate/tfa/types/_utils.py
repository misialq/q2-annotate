# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from urllib.parse import quote


def _gene_taxonomy_id(taxon_id: str, gene_id: str) -> str:
    """Join readable IDs, escaping the separator and literal percent signs."""
    return f"{quote(str(taxon_id), safe='')}|{quote(str(gene_id), safe='')}"
