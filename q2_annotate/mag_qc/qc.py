# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import shutil
from importlib import resources

import pandas as pd
import q2templates
from q2_types.kraken2 import Kraken2OutputDirectoryFormat
import argparse
import json
import os
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from itertools import takewhile
from typing import List, Dict
import sys
import re

import skbio.io
import multiprocessing
from skbio import DNA

from q2_types.per_sample_sequences import MultiFASTADirectoryFormat

from q2_annotate.kraken2.select import _find_lcas


TEMPLATES = resources.files("q2_annotate") / "assets"


def parse_busco_file(df: pd.DataFrame) -> dict:
    """Parse BUSCO results file to extract MAG quality metrics."""
    busco_data = {}
    for mag_id, row in df.iterrows():
        busco_data[mag_id] = {
            'sample_id': row['sample_id'],
            'completeness': row['completeness'],
            'contamination': row['contamination']
        }
    return busco_data


def parse_taxonomy_file(taxonomy: pd.Series):
    """Parse taxonomy mapping file to map taxon IDs to full taxonomy."""
    # Create mapping from taxon ID to taxonomy
    taxonomy_map = {}
    for taxon_id, taxon in taxonomy.iteritems():
        taxonomy_map[taxon_id] = taxon

    print(f"Loaded {len(taxonomy_map)} taxonomy mappings")
    return taxonomy_map


def extract_assigned_taxonomy_lca(mag_id, contigs_data):
    """Extract the assigned taxonomy for a MAG using LCA calculation."""
    if not contigs_data:
        return "Unassigned"

    # Prepare data for LCA calculation
    taxa_data = []

    for taxonomy in contigs_data.keys():
        for contig in contigs_data[taxonomy]:
            taxa_data.append({
                'mag_id': mag_id,
                'Taxon': taxonomy
            })

    if not taxa_data:
        return "Unassigned"

    # Convert to DataFrame
    taxa_df = pd.DataFrame(taxa_data)

    # Calculate LCA
    try:
        lca_result = _find_lcas([taxa_df], mode='lca')
        if mag_id in lca_result.index:
            lca_taxonomy = lca_result.loc[mag_id, 'Taxon']
            return lca_taxonomy if isinstance(lca_taxonomy, str) else "Unassigned"

    except Exception as e:
        print(f"Warning: LCA calculation failed for MAG {mag_id}: {e}")
        # Fallback to most abundant genus
        return extract_assigned_taxonomy_fallback(contigs_data)

    return "Unassigned"


def extract_assigned_taxonomy_fallback(contigs_data):
    """Fallback method: Extract the assigned taxonomy based on most abundant genus."""
    if not contigs_data:
        return "Unassigned"

    # Count taxa at genus level (or highest available level)
    genus_counts = defaultdict(int)

    for taxonomy in contigs_data.keys():
        # Split taxonomy and find genus level
        levels = taxonomy.split(';')
        genus_level = None
        for level in levels:
            if level.startswith('g__'):
                genus_level = level
                break

        if genus_level:
            genus_counts[genus_level] += len(contigs_data[taxonomy])

    if genus_counts:
        # Return the most abundant genus
        most_abundant_genus = max(genus_counts, key=genus_counts.get)
        return most_abundant_genus

    return "Unassigned"


def parse_coverage_file(df: pd.DataFrame):
    """Parse coverage TSV file to get coverage data per contig per sample."""
    coverage_data = {}

    for sample_id, row in df.iterrows():
        for contig_id in df.columns:
            value = row[contig_id]
            if value is not None and value > 0:
                coverage_data[contig_id] = value
    return coverage_data


def parse_kraken2_output(output_file, taxonomy_map):
    """Parse a single Kraken2 output file to get taxonomic assignments."""
    contigs_data = defaultdict(list)

    with open(output_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 4:
                continue

            classification = parts[0]  # C = classified, U = unclassified
            contig_id = parts[1]  # Contig/sequence ID
            taxon_id = parts[2]  # Taxon ID

            # Only process classified sequences
            if classification == 'C' and taxon_id in taxonomy_map:
                taxonomy = taxonomy_map[taxon_id]

                # Group by taxonomy at reasonable level (up to genus/species)
                contigs_data[taxonomy].append({
                    'contig_id': contig_id,
                    'specific_taxonomy': taxonomy,
                })

    return dict(contigs_data)


def _parse_single_kraken(args):
    mag_id, sample_id, fp, taxonomy_map = args
    taxonomy = parse_kraken2_output(fp, taxonomy_map)
    return mag_id, sample_id, taxonomy


def parse_kraken2_directory(kraken2_dir: Kraken2OutputDirectoryFormat, taxonomy_map: dict, cpus: int):
    """Parse directory structure containing Kraken2 output files."""
    mag_taxonomies = {}
    tasks = []
    for sample_id, mag_dict in kraken2_dir.file_dict().items():
        for mag_id, fp in mag_dict.items():
            tasks.append((mag_id, sample_id, fp, taxonomy_map))

    with multiprocessing.Pool(processes=cpus) as pool:
        for mag_id, sample_id, taxonomy in pool.map(_parse_single_kraken, tasks):
            mag_taxonomies[mag_id] = {
                'sample_id': sample_id,
                'taxonomy': taxonomy
            }

    return mag_taxonomies


def calculate_mag_stats(contigs_data):
    """Calculate MAG statistics from contig data."""
    if not contigs_data:
        return {
            'n50': 0,
            'l50': 0,
            'total_length': 0,
            'quality': 'Low'
        }

    # Collect all contig lengths
    all_lengths = []
    total_length = 0

    for taxonomy_group in contigs_data.values():
        for contig in taxonomy_group:
            length = contig['length']
            all_lengths.append(length)
            total_length += length

    # Sort lengths in descending order for N50 calculation
    all_lengths.sort(reverse=True)

    # Calculate N50 and L50
    n50 = 0
    l50 = 0
    cumulative_length = 0
    half_total = total_length / 2

    for i, length in enumerate(all_lengths):
        cumulative_length += length
        if cumulative_length >= half_total:
            n50 = length
            l50 = i + 1
            break

    return {
        'n50': n50,
        'l50': l50,
        'total_length': total_length
    }


def determine_quality(busco_results: pd.DataFrame):
    """Determine MAG quality based on completeness and contamination."""
    completeness = float(busco_results['completeness'])
    contamination = float(busco_results['contamination'])
    if completeness >= 90 and contamination <= 5:
        return "High"
    elif completeness >= 50 and contamination <= 10:
        return "Medium"
    else:
        return "Low"


def clean_nan_values(obj):
    """Recursively clean NaN values from nested data structures, replacing with None."""
    if isinstance(obj, dict):
        return {key: clean_nan_values(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [clean_nan_values(item) for item in obj]
    elif pd.isna(obj):
        return None
    elif isinstance(obj, (np.integer, np.floating)):
        if pd.isna(obj):
            return None
        return obj.item()  # Convert numpy types to Python native types
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj


def add_contig_metrics(
        contigs_data: dict, contig_metrics: dict
) -> dict:
    """Add GC content, length and coverage for each contig of every MAG.

        Args:
            contigs_data (dict): A dictionary containing taxonoym data per MAG.
            mags (MultiFASTADirectoryFormat): MAG fasta files.
            coverage: A DataFrame of contig coverage.

        Returns:
            A dictionary where keys are MAG IDs (derived from file names)
            and values are dictionaries mapping contig IDs to their GC content.
            Example: {mag_id: {contig_id: gc_content, coverage: depth}, ...}
        """
    results = {}
    for taxonomy, contig_list in contigs_data.items():
        contigs = []
        for contig in contig_list:
            _id = contig['contig_id']
            contig.update(contig_metrics[_id])
            contigs.append(contig)
        results[taxonomy] = contigs
    return results



def _process_single_fasta(args):
    fp, coverage = args
    file_results = {}
    for seq in skbio.io.read(fp, format='fasta'):
        contig_id = seq.metadata['id']
        seq = DNA(seq)
        file_results[contig_id] = {
            "gc_content": round(100 * seq.gc_content(), 2),
            "length": len(seq),
            "coverage": round(coverage.get(contig_id, 0), 2)
        }
    return file_results

def calculate_contig_metrics(mags: MultiFASTADirectoryFormat, coverage: dict, cpus: int) -> dict:
    """Add GC content, length and coverage for each contig of every MAG.

        Args:
            mags (MultiFASTADirectoryFormat): MAG fasta files.

        Returns:
            A dictionary where keys are MAG IDs (derived from file names)
            and values are dictionaries mapping contig IDs to their GC content.
            Example: {mag_id: {contig_id: gc_content, coverage: depth}, ...}
        """
    results = {}
    tasks = []
    for sample_id, mag_dict in mags.sample_dict().items():
        for _, fp in mag_dict.items():
            tasks.append((fp, coverage))

    with multiprocessing.Pool(processes=cpus) as pool:
        for file_result in pool.map(_process_single_fasta, tasks):
            results.update(file_result)

    return results


def visualize_mag_quality(
        output_dir: str,
        mags: MultiFASTADirectoryFormat,
        busco_results: pd.DataFrame,
        kraken2_outputs: Kraken2OutputDirectoryFormat,
        mag_taxonomy: pd.Series,
        contig_coverage: pd.DataFrame,
        n_cpus: int = 1
):
    # Parse input files
    busco_data = parse_busco_file(busco_results)
    taxonomy_map = mag_taxonomy.to_dict()
    coverage_data = parse_coverage_file(contig_coverage)
    mag_taxonomies = parse_kraken2_directory(kraken2_outputs, taxonomy_map, n_cpus)
    contig_metrics = calculate_contig_metrics(mags, coverage_data, n_cpus)

    result = defaultdict(list)
    for mag_id, busco_info in busco_data.items():
        sample_id = busco_info['sample_id']

        # Get taxonomic data for this MAG
        taxonomy_info = mag_taxonomies.get(mag_id, {
            'sample_id': sample_id,
            'taxonomy': {}
        })

        # Add contig length, GC content and coverage
        contigs_data = add_contig_metrics(taxonomy_info['taxonomy'], contig_metrics)

        # Determine assigned taxonomy using LCA
        assigned_taxonomy = extract_assigned_taxonomy_lca(mag_id, contigs_data)

        # Determine quality
        quality = determine_quality(busco_info)

        # Create MAG entry
        mag_entry = {
            'id': mag_id,
            'contamination': busco_info['contamination'],
            'completeness': busco_info['completeness'],
            'assigned_taxonomy': assigned_taxonomy,
            'contigs': contigs_data,
            'quality': quality,
            **calculate_mag_stats(contigs_data)
        }

        result[sample_id].append(mag_entry)

    # Clean NaN values before JSON serialization
    cleaned_result = clean_nan_values(dict(result))

    with open(os.path.join(output_dir, "data.json"), 'w') as f:
        json.dump(cleaned_result, f, indent=2)

    templates = [TEMPLATES / "mag_qc" / "index.html"]

    # Copy JS/CSS files
    for fn in ("scripts.js", "bootstrapMagic.js", "styles.css"):
        shutil.copy(str(TEMPLATES / "mag_qc" / fn), os.path.join(output_dir, fn))

    q2templates.render(templates, output_dir, context={})
