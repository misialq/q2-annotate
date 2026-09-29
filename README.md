# q2-annotate
![CI](https://github.com/bokulich-lab/q2-annotate/actions/workflows/ci.yaml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-annotate/graph/badge.svg?token=PSCAYJUP01)](https://codecov.io/gh/bokulich-lab/q2-annotate)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QIIME 2 plugin for functional annotation and taxonomic classification of shotgun metagenomes.

## Taxonomic functional attribution

`estimate-tfa` returns two artifacts:

- `feature_load`: a sparse `FeatureTable[Frequency]` with the original sample
  IDs as columns and one feature per observed taxon/gene pair.
- `gene_taxonomy`: `FeatureData[GeneTaxonomy]`, a readable TSV mapping each
  `Feature ID` to its `Taxon ID`, `Gene ID`, and `Taxon` label.

Feature IDs are `gt_` followed by the full SHA-256 digest of the UTF-8 JSON
array `[taxon ID, gene ID]` (compact separators, no ASCII escaping). They are
stable across runs and input ordering, and do not change when other pairs are
added. The mapping preserves the original IDs without delimiter ambiguity.

For taxon *t*, gene *g*, and sample *s*, a load is the sum over contigs assigned
to *t* of contig abundance in *s* multiplied by the contig's count of *g*.
Only pairs present in the gene inventory become rows; sample zeros stay
sparse. Fractional abundance-weighted values are preserved without rounding,
following the `FeatureTable[Frequency]` convention used by `q2-mag
estimate-abundance` for RPKM/TPM. The obsolete `TFA` semantic type and its
custom BIOM format have been removed; both outputs must be retained to
interpret the table's feature IDs.

The sample-resolved output replaces the older `estimate-tfa` layout, which
summed abundances across samples and used functional feature IDs as sample IDs.
Run `estimate-tfa` again from the original inputs to obtain per-sample loads;
the original layout cannot recover them.

## Installation
_q2-annotate_ is available as part of the QIIME 2 moshpit distribution. For installation and usage instructions please consult the official [QIIME 2 documentation](https://docs.qiime2.org).

## Functionality
This QIIME 2 plugin contains actions used to annotate and classify (meta)genomes:

| Action                    | Description                                                                            | Underlying tool       |
|---------------------------|----------------------------------------------------------------------------------------|-----------------------|
| build-custom-diamond-db   | Create a DIAMOND reference database from a FASTA input file.                           | [Diamond](https://github.com/bbuchfink/diamond) |
| build-eggnog-diamond-db   | Create a DIAMOND reference database for the specified taxon.                           | [Diamond](https://github.com/bbuchfink/diamond) |
| build-kraken-db           | Fetch an existing or build a custom Kraken 2 database.                                 | [Kraken 2](https://ccb.jhu.edu/software/kraken2/)      |
| classify-kaiju            | Classify reads using Kaiju.                                                            | [Kaiju](https://bioinformatics-centre.github.io/kaiju/) |
| classify-kraken2          | Classify reads/MAGs using Kraken 2.                                                    | [Kraken 2](https://ccb.jhu.edu/software/kraken2/)      |
| construct-pangenome-index | Construct the human pangenome Bowtie 2 index.                                          | [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) |
| estimate-bracken          | Perform read abundance re-estimation using Bracken.                                    | [Kraken 2](https://ccb.jhu.edu/software/bracken/) |
| estimate-tfa              | Estimate a sparse load table and gene taxonomy mapping. | - |
| extract-annotations       | Extract annotation frequencies from all annotations.                                   | - |
| fetch-diamond-db          | Fetch the complete Diamond database necessary to run the eggnog-diamond-search action. | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-eggnog-db           | Fetch the databases necessary to run the eggnog-annotate action.                       | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-eggnog-hmmer-db     | Fetch the taxon specific database necessary to run the eggnog-hmmer-search action.     | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-eggnog-proteins     | Fetch the databases necessary to run the build-eggnog-diamond-db action.               | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| fetch-kaiju-db            | Fetch Kaiju database.                                                                  | [Kaiju](https://bioinformatics-centre.github.io/kaiju/) |
| fetch-ncbi-taxonomy       | Fetch NCBI reference taxonomy.                                                         | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| filter-kraken2-reads      | Filter Kraken 2-classified reads by taxonomy (for example, retain only bacterial reads). | - |
| filter-kraken2-results    | Filter Kraken 2 reports and outputs by metadata and/or abundance                       | - |
| filter-reads-pangenome    | Remove contaminating human reads.                                                      | [Bowtie 2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) |
| inspect-kraken2-db        | Inspect a Kraken 2 database.                                                           | - |
| kraken2-to-features       | Select downstream features from Kraken 2.                                              | - |
| kraken2-to-mag-features   | Select downstream MAG features from Kraken 2.                                          | - |
| map-eggnog                | Annotate orthologs against eggNOG database.                                            | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| multiply-tables           | Multiply two feature tables.                                                           | - |
| predict-genes-prodigal    | Predict gene sequences from MAGs using Prodigal.                                       | [Prodigal](https://github.com/hyattpd/Prodigal) |
| search-orthologs-diamond  | Run eggNOG search using Diamond aligner.                                               | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
| search-orthologs-hmmer    | Run eggNOG search using HMMER aligner.                                                 | [EggNOG mapper](https://github.com/eggnogdb/eggnog-mapper) |
