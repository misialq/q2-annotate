# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import csv

from qiime2.plugin import ValidationError, model

GENE_TAXONOMY_COLUMNS = ["Taxon ID", "Gene ID", "Taxon"]


class GeneTaxonomyFormat(model.TextFileFormat):
    """TSV mapping ordinary feature IDs to taxon IDs, gene IDs, and labels."""

    def _validate_(self, level):
        seen_ids, seen_pairs = set(), set()
        try:
            with self.open() as handle:
                reader = csv.reader(handle, delimiter="\t", strict=True)
                if next(reader, None) != ["Feature ID", *GENE_TAXONOMY_COLUMNS]:
                    raise ValidationError(
                        "Expected Feature ID, Taxon ID, Gene ID, Taxon header."
                    )
                # Minimal validation checks the first 100 records; full validation streams all.
                for index, row in enumerate(reader):
                    if level == "min" and index >= 100:
                        break
                    if len(row) != 4 or any(not value.strip() for value in row[:3]):
                        raise ValidationError(
                            "Each gene taxonomy row requires four fields and nonempty IDs."
                        )
                    pair = tuple(row[1:3])
                    # Both the feature ID and taxon/gene pair must identify one row.
                    if row[0] in seen_ids or pair in seen_pairs:
                        raise ValidationError(
                            "Gene taxonomy feature IDs and taxon/gene pairs must be unique."
                        )
                    seen_ids.add(row[0])
                    seen_pairs.add(pair)
        except (UnicodeDecodeError, csv.Error) as error:
            raise ValidationError("Expected a UTF-8 gene taxonomy TSV file.") from error


GeneTaxonomyDirFmt = model.SingleFileDirectoryFormat(
    "GeneTaxonomyDirFmt", "gene-taxonomy.tsv", GeneTaxonomyFormat
)
