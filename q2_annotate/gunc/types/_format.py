# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin import model


class GUNCResultsFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


class GUNCResultsDirectoryFormat(model.SingleFileDirectoryFormat):
    def __init__(self, path=None, mode="w"):
        super().__init__(
            path, mode,
            file_name="gunc_results.tsv",
            format=GUNCResultsFormat,
        )


class GUNCDatabaseDirFmt(model.DirectoryFormat):
    def _validate_(self, level):
        pass
