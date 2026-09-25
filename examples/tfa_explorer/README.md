# Synthetic TFA explorer dataset

This example is **simulated**, not derived from biological samples. It is
designed to exercise the TFA visualizer with 60 samples (30 paired subjects,
baseline and day 14), 13 taxon IDs, 14 functions, and 58 observed
taxon–function pairs. Fifteen subjects are assigned to each of two synthetic
treatment groups. The generator introduces a larger day 14 resistance signal
in selected taxa in the antibiotic group, along with lower simulated loads
for some commensal taxa. These patterns are illustration only and should not
be interpreted as biological findings.

The source files are readable TSVs:

- `tfa-table.tsv`: one row per taxon–function pair, one column per sample;
  zeros represent absent loads. The corresponding QIIME artifact is
  `feature-load.qza` (`FeatureTable[TFA]`).
- `taxonomy.tsv`: `FeatureData[Taxonomy]` assignments for 12 of the 13 taxon
  IDs. `T013` is intentionally unassigned to demonstrate ID fallback. The
  corresponding artifact is `taxonomy.qza`.
- `sample-metadata.tsv`: subject, treatment, visit, combined treatment/visit,
  batch, and age. The categorical columns appear in the visualizer's group
  selector.

With this branch installed in a QIIME 2 environment, run from the repository
root:

```bash
qiime annotate explore-tfa \
  --i-feature-load examples/tfa_explorer/feature-load.qza \
  --i-taxonomy examples/tfa_explorer/taxonomy.qza \
  --m-metadata-file examples/tfa_explorer/sample-metadata.tsv \
  --o-visualization tfa-explorer-demo.qzv
```

Open `tfa-explorer-demo.qzv` with `qiime tools view`. Select
`treatment_visit` to compare the four groups. Try `blaCTX-M` and
`butyryl_CoA_transferase` in the function selector, then select individual
taxa or click heatmap cells. Omit `--i-taxonomy` to see taxon IDs.

To regenerate all TSV and QZA files from the fixed random seed, run
`python examples/tfa_explorer/generate_demo.py` in the same environment.
