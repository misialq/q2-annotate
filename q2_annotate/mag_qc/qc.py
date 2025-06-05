# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import shutil
from importlib import resources
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple
import json
import os
import multiprocessing

import pandas as pd
import numpy as np
import q2templates
import skbio.io
from skbio import DNA
from collections import defaultdict

from q2_types.kraken2 import Kraken2OutputDirectoryFormat
from q2_types.per_sample_sequences import MultiFASTADirectoryFormat
from q2_annotate.kraken2.select import _find_lcas

TEMPLATES = resources.files("q2_annotate") / "assets"


@dataclass
class ContigData:
    """Single contig information."""
    contig_id: str
    length: int
    gc_content: float
    coverage: float
    taxonomy: str = "Unassigned"


@dataclass
class MAGData:
    """Complete MAG information."""
    id: str
    sample_id: str
    completeness: float
    contamination: float
    quality: str
    assigned_taxonomy: str
    contigs: Dict[str, List[dict]]
    n50: int = 0
    l50: int = 0
    total_length: int = 0


def determine_quality(completeness: float, contamination: float) -> str:
    """Determine MAG quality based on completeness and contamination."""
    if completeness >= 90 and contamination <= 5:
        return "High"
    elif completeness >= 50 and contamination <= 10:
        return "Medium"
    else:
        return "Low"


def calculate_assembly_stats(contigs_by_taxonomy: Dict[str, List[dict]]) -> Tuple[int, int, int]:
    """Calculate N50, L50, and total length from contig data."""
    if not contigs_by_taxonomy:
        return 0, 0, 0
    
    # Collect all contig lengths from all taxonomies
    lengths = []
    for taxonomy_contigs in contigs_by_taxonomy.values():
        lengths.extend([contig['length'] for contig in taxonomy_contigs])
    
    if not lengths:
        return 0, 0, 0
    
    lengths.sort(reverse=True)
    total_length = sum(lengths)
    
    n50 = l50 = 0
    cumulative = 0
    half_total = total_length / 2
    
    for i, length in enumerate(lengths):
        cumulative += length
        if cumulative >= half_total:
            n50 = length
            l50 = i + 1
            break
    
    return n50, l50, total_length


def extract_assigned_taxonomy_lca(contigs_by_taxonomy: Dict[str, List[dict]]) -> str:
    """Extract assigned taxonomy using LCA calculation."""
    if not contigs_by_taxonomy:
        return "Unassigned"
    
    # Collect all assigned taxonomies (excluding "Unassigned")
    assigned_taxonomies = [taxonomy for taxonomy in contigs_by_taxonomy.keys() 
                          if taxonomy != "Unassigned"]
    
    if not assigned_taxonomies:
        return "Unassigned"
    
    # Prepare data for LCA calculation - weight by number of contigs
    taxa_data = []
    for taxonomy in assigned_taxonomies:
        for _ in contigs_by_taxonomy[taxonomy]:  # Add one entry per contig
            taxa_data.append({'mag_id': 'temp', 'Taxon': taxonomy})
    
    if not taxa_data:
        return "Unassigned"
    
    taxa_df = pd.DataFrame(taxa_data)
    
    try:
        lca_result = _find_lcas([taxa_df], mode='lca')
        if 'temp' in lca_result.index:
            return lca_result.loc['temp', 'Taxon']
    except Exception as e:
        print(f"Warning: LCA calculation failed: {e}")
        # Fallback to most common taxonomy (by number of contigs)
        taxonomy_counts = {taxonomy: len(contigs) 
                          for taxonomy, contigs in contigs_by_taxonomy.items()
                          if taxonomy != "Unassigned"}
        if taxonomy_counts:
            return max(taxonomy_counts, key=taxonomy_counts.get)
    
    return "Unassigned"


def _process_single_mag(args) -> Tuple[str, str, Dict[str, List[dict]]]:
    """Process a single MAG file to extract contig metrics and taxonomy."""
    mag_id, sample_id, mag_fp, kraken_fp, coverage_data, taxonomy_map = args
    
    contigs_by_taxonomy = defaultdict(list)
    
    # Load taxonomy assignments for this MAG
    contig_taxonomy = {}
    if kraken_fp and os.path.exists(kraken_fp):
        with open(kraken_fp, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3 and parts[0] == 'C':  # Classified
                    contig_id, taxon_id = parts[1], parts[2]
                    if taxon_id in taxonomy_map:
                        contig_taxonomy[contig_id] = taxonomy_map[taxon_id]
    
    # Process FASTA file
    for seq in skbio.io.read(mag_fp, format='fasta'):
        contig_id = seq.metadata['id']
        seq_dna = DNA(seq)
        
        taxonomy = contig_taxonomy.get(contig_id, "Unassigned")
        
        contig_data = {
            'contig_id': contig_id,
            'length': len(seq_dna),
            'gc_content': round(100 * seq_dna.gc_content(), 2),
            'coverage': round(coverage_data.get(contig_id, 0), 2),
            'taxonomy': taxonomy
        }
        
        contigs_by_taxonomy[taxonomy].append(contig_data)
    
    return mag_id, sample_id, dict(contigs_by_taxonomy)


def clean_for_json(obj):
    """Clean data for JSON serialization."""
    if isinstance(obj, dict):
        return {k: clean_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [clean_for_json(item) for item in obj]
    elif pd.isna(obj):
        return None
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def visualize_mag_quality(
    output_dir: str,
    mags: MultiFASTADirectoryFormat,
    busco_results: pd.DataFrame,
    kraken2_outputs: Kraken2OutputDirectoryFormat,
    mag_taxonomy: pd.Series,
    contig_coverage: pd.DataFrame,
    n_cpus: int = 1
):
    """Simplified main function with streamlined data processing."""
    
    # Prepare data
    taxonomy_map = mag_taxonomy.to_dict()
    
    # Parse coverage data - Fix the original parsing logic
    coverage_data = {}
    for contig_id in contig_coverage.columns:
        for sample_id in contig_coverage.index:
            value = contig_coverage.loc[sample_id, contig_id]
            if pd.notna(value) and value > 0:
                coverage_data[contig_id] = float(value)
    
    # Prepare tasks for parallel processing
    tasks = []
    kraken_files = kraken2_outputs.file_dict()
    
    for sample_id, mag_dict in mags.sample_dict().items():
        for mag_id, mag_fp in mag_dict.items():
            # Find corresponding Kraken2 file
            kraken_fp = None
            if sample_id in kraken_files and mag_id in kraken_files[sample_id]:
                kraken_fp = kraken_files[sample_id][mag_id]
            
            tasks.append((mag_id, sample_id, mag_fp, kraken_fp, coverage_data, taxonomy_map))
    
    # Process all MAGs in parallel
    results = {}
    with multiprocessing.Pool(processes=n_cpus) as pool:
        for mag_id, sample_id, contigs_by_taxonomy in pool.map(_process_single_mag, tasks):
            if sample_id not in results:
                results[sample_id] = []
            
            # Get BUSCO data for this MAG
            if mag_id in busco_results.index:
                completeness = float(busco_results.loc[mag_id, 'completeness'])
                contamination = float(busco_results.loc[mag_id, 'contamination'])
                quality = determine_quality(completeness, contamination)
                
                # Calculate assembly stats
                n50, l50, total_length = calculate_assembly_stats(contigs_by_taxonomy)
                
                # Determine assigned taxonomy
                assigned_taxonomy = extract_assigned_taxonomy_lca(contigs_by_taxonomy)
                
                # Create MAG object
                mag = MAGData(
                    id=mag_id,
                    sample_id=sample_id,
                    completeness=completeness,
                    contamination=contamination,
                    quality=quality,
                    assigned_taxonomy=assigned_taxonomy,
                    contigs=contigs_by_taxonomy,
                    n50=n50,
                    l50=l50,
                    total_length=total_length
                )
                
                results[sample_id].append(asdict(mag))
    
    # Clean and save results
    cleaned_results = clean_for_json(results)
    
    with open(os.path.join(output_dir, "data.json"), 'w') as f:
        json.dump(cleaned_results, f, indent=2)
    
    # Copy template files
    templates = [TEMPLATES / "mag_qc" / "index.html"]
    for fn in ("scripts.js", "bootstrapMagic.js", "styles.css"):
        shutil.copy(str(TEMPLATES / "mag_qc" / fn), os.path.join(output_dir, fn))
    
    q2templates.render(templates, output_dir, context={})
