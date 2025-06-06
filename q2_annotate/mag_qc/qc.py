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
import logging
import time

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

# Configure logging for performance monitoring
def configure_logging(level: str = "INFO", 
                     format_string: str = None,
                     enable_timing: bool = True) -> None:
    """Configure logging for MAG QC performance monitoring.
    
    Args:
        level: Logging level ('DEBUG', 'INFO', 'WARNING', 'ERROR')
        format_string: Custom format string for log messages
        enable_timing: Whether to include timing information in logs
    """
    if format_string is None:
        if enable_timing:
            format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        else:
            format_string = "%(name)s - %(levelname)s - %(message)s"
    
    # Configure the logger for this module
    logger = logging.getLogger(__name__)
    
    # Remove any existing handlers to avoid duplicates
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Create console handler
    handler = logging.StreamHandler()
    handler.setLevel(getattr(logging, level.upper()))
    
    # Create formatter
    formatter = logging.Formatter(format_string)
    handler.setFormatter(formatter)
    
    # Add handler to logger
    logger.addHandler(handler)
    logger.setLevel(getattr(logging, level.upper()))
    
    # Prevent propagation to avoid duplicate messages
    logger.propagate = False

# Set up default logging configuration
configure_logging(level="DEBUG", enable_timing=True)

# Get the module logger
logger = logging.getLogger(__name__)


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
    start_time = time.time()
    logger.debug("Starting assembly stats calculation")
    
    if not contigs_by_taxonomy:
        logger.debug("No contigs provided for assembly stats")
        return 0, 0, 0
    
    # Use generator to avoid creating intermediate lists
    logger.debug("Collecting contig lengths from all taxonomies")
    lengths = [contig['length'] 
              for taxonomy_contigs in contigs_by_taxonomy.values()
              for contig in taxonomy_contigs]
    
    if not lengths:
        logger.debug("No contig lengths found")
        return 0, 0, 0
    
    logger.debug(f"Processing {len(lengths)} contigs for assembly stats")
    
    # Sort and calculate total in one pass
    sort_start = time.time()
    lengths.sort(reverse=True)
    logger.debug(f"Sorting completed in {time.time() - sort_start:.3f} seconds")
    
    total_length = sum(lengths)
    half_total = total_length / 2
    logger.debug(f"Total assembly length: {total_length:,} bp")
    
    # Calculate N50/L50 in single loop
    cumulative = 0
    for i, length in enumerate(lengths):
        cumulative += length
        if cumulative >= half_total:
            n50, l50 = length, i + 1
            total_time = time.time() - start_time
            logger.debug(f"Assembly stats calculated in {total_time:.3f} seconds: N50={n50:,}, L50={l50}, Total={total_length:,}")
            return n50, l50, total_length
    
    # Fallback (shouldn't reach here with valid data)
    logger.warning("Assembly stats calculation reached fallback case")
    return 0, 0, total_length


def extract_assigned_taxonomy_lca(contigs_by_taxonomy: Dict[str, List[dict]]) -> str:
    """Extract assigned taxonomy using LCA calculation."""
    start_time = time.time()
    logger.debug("Starting LCA taxonomy extraction")
    
    if not contigs_by_taxonomy:
        logger.debug("No contigs provided for LCA calculation")
        return "Unassigned"
    
    # Collect all assigned taxonomies (excluding "Unassigned")
    assigned_taxonomies = [taxonomy for taxonomy in contigs_by_taxonomy.keys() 
                          if taxonomy != "Unassigned"]
    
    if not assigned_taxonomies:
        logger.debug("No assigned taxonomies found, returning Unassigned")
        return "Unassigned"
    
    logger.debug(f"Found {len(assigned_taxonomies)} assigned taxonomies for LCA calculation")
    
    # Prepare data for LCA calculation - weight by number of contigs
    logger.debug("Preparing taxa data for LCA calculation")
    taxa_prep_start = time.time()
    taxa_data = [{'mag_id': 'temp', 'Taxon': taxonomy}
                 for taxonomy in assigned_taxonomies
                 for _ in contigs_by_taxonomy[taxonomy]]
    logger.debug(f"Taxa data preparation completed in {time.time() - taxa_prep_start:.3f} seconds, {len(taxa_data)} entries")
    
    if not taxa_data:
        logger.debug("No taxa data prepared, returning Unassigned")
        return "Unassigned"
    
    df_start = time.time()
    taxa_df = pd.DataFrame(taxa_data)
    logger.debug(f"DataFrame creation completed in {time.time() - df_start:.3f} seconds")
    
    try:
        logger.debug("Executing LCA calculation")
        lca_start = time.time()
        lca_result = _find_lcas([taxa_df], mode='lca')
        logger.debug(f"LCA calculation completed in {time.time() - lca_start:.3f} seconds")
        
        if 'temp' in lca_result.index:
            result_taxonomy = lca_result.loc['temp', 'Taxon']
            total_time = time.time() - start_time
            logger.info(f"LCA taxonomy extraction completed in {total_time:.3f} seconds: {result_taxonomy}")
            return result_taxonomy
    except Exception as e:
        logger.warning(f"LCA calculation failed: {e}, falling back to most common taxonomy")
        # Fallback to most common taxonomy (by number of contigs)
        taxonomy_counts = {taxonomy: len(contigs) 
                          for taxonomy, contigs in contigs_by_taxonomy.items()
                          if taxonomy != "Unassigned"}
        if taxonomy_counts:
            fallback_result = max(taxonomy_counts, key=taxonomy_counts.get)
            total_time = time.time() - start_time
            logger.info(f"LCA fallback completed in {total_time:.3f} seconds: {fallback_result}")
            return fallback_result
    
    logger.debug("LCA calculation returned no result, returning Unassigned")
    return "Unassigned"


def _process_single_mag(args) -> Tuple[str, str, Dict[str, List[dict]]]:
    """Process a single MAG file to extract contig metrics and taxonomy."""
    start_time = time.time()
    mag_id, sample_id, mag_fp, kraken_fp, coverage_data, taxonomy_map = args
    
    logger.debug(f"Starting processing of MAG {mag_id} from sample {sample_id}")
    
    contigs_by_taxonomy = defaultdict(list)
    
    # Load taxonomy assignments for this MAG
    contig_taxonomy = {}
    if kraken_fp and os.path.exists(kraken_fp):
        logger.debug(f"Loading taxonomy assignments from {kraken_fp}")
        kraken_start = time.time()
        line_count = 0
        classified_count = 0
        
        with open(kraken_fp, 'r') as f:
            for line in f:
                line_count += 1
                if line.startswith('C\t'):  # Fast prefix check before splitting
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 3:
                        contig_id, taxon_id = parts[1], parts[2]
                        # Direct lookup without 'in' check (KeyError handling below)
                        try:
                            contig_taxonomy[contig_id] = taxonomy_map[taxon_id]
                            classified_count += 1
                        except KeyError:
                            pass  # Taxon ID not in map, skip
        
        kraken_time = time.time() - kraken_start
        logger.debug(f"Kraken file processing completed in {kraken_time:.3f} seconds: {line_count} lines, {classified_count} classified contigs")
    else:
        logger.debug(f"No Kraken file found for MAG {mag_id}, all contigs will be unassigned")
    
    # Process FASTA file
    logger.debug(f"Processing FASTA file {mag_fp}")
    fasta_start = time.time()
    contig_count = 0
    
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
        contig_count += 1
    
    fasta_time = time.time() - fasta_start
    total_time = time.time() - start_time
    
    logger.debug(f"FASTA processing completed in {fasta_time:.3f} seconds: {contig_count} contigs")
    logger.info(f"MAG {mag_id} processing completed in {total_time:.3f} seconds: {contig_count} contigs, {len(contigs_by_taxonomy)} taxonomies")
    
    return mag_id, sample_id, dict(contigs_by_taxonomy)


def clean_for_json(obj):
    """Clean data for JSON serialization with optimized type checking."""
    # Check numpy types first (most common in scientific data)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif isinstance(obj, dict):
        return {k: clean_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [clean_for_json(item) for item in obj]
    # Check for scalar NaN values (avoid calling pd.isna on collections)
    elif hasattr(obj, '__len__') and not isinstance(obj, (str, bytes)):
        # This is a collection type, don't check pd.isna
        return obj
    elif pd.isna(obj):
        return None
    return obj


def parse_coverage_data(contig_coverage: pd.DataFrame) -> dict:
    """Parse coverage DataFrame to get coverage data per contig.
    
    Assumes each contig appears in only one sample (i.e., has positive coverage 
    in only one row). This assumption enables highly optimized processing using 
    vectorized pandas operations.
    
    Args:
        contig_coverage: DataFrame with samples as rows and contigs as columns
        
    Returns:
        Dictionary mapping contig IDs to coverage values
    """
    start_time = time.time()
    logger.info(f"Starting coverage data parsing for {contig_coverage.shape[0]} samples and {contig_coverage.shape[1]} contigs")
    
    # Early exit for empty DataFrame
    if contig_coverage.empty:
        logger.warning("Empty coverage DataFrame provided")
        return {}
    
    # Optimized path: each contig is present in only one sample
    logger.debug("Converting DataFrame to stacked Series")
    stack_start = time.time()
    stacked = contig_coverage.stack()
    logger.debug(f"DataFrame stacking completed in {time.time() - stack_start:.3f} seconds")
    
    logger.debug("Filtering positive coverage values")
    filter_start = time.time()
    positive_values = stacked[stacked > 0]
    logger.debug(f"Filtering completed in {time.time() - filter_start:.3f} seconds, found {len(positive_values)} positive values")
    
    # Convert directly to dict (contig_id is now the index level 1)
    logger.debug("Converting to dictionary")
    dict_start = time.time()
    coverage_data = {contig_id: float(value) 
                    for (sample_id, contig_id), value in positive_values.items()}
    logger.debug(f"Dictionary conversion completed in {time.time() - dict_start:.3f} seconds")
    
    total_time = time.time() - start_time
    logger.info(f"Coverage data parsing completed in {total_time:.3f} seconds, extracted {len(coverage_data)} contigs")
    return coverage_data


def process_mags_parallel(
    mags: MultiFASTADirectoryFormat,
    busco_results: pd.DataFrame,
    kraken2_outputs: Kraken2OutputDirectoryFormat,
    taxonomy_map: dict,
    coverage_data: dict,
    n_cpus: int = 1
) -> dict:
    """Process all MAGs in parallel and return organized results.
    
    Args:
        mags: MAG sequences directory
        busco_results: DataFrame with MAG quality metrics
        kraken2_outputs: Kraken2 taxonomic classification results
        taxonomy_map: Mapping from taxon IDs to taxonomy strings
        coverage_data: Contig coverage information
        n_cpus: Number of CPU cores to use
        
    Returns:
        Dictionary with sample_id as keys and lists of MAG data as values
    """
    start_time = time.time()
    logger.info(f"Starting parallel MAG processing with {n_cpus} CPUs")
    
    # Prepare tasks for parallel processing
    logger.debug("Preparing tasks for parallel processing")
    task_prep_start = time.time()
    tasks = []
    kraken_files = kraken2_outputs.file_dict()
    
    for sample_id, mag_dict in mags.sample_dict().items():
        for mag_id, mag_fp in mag_dict.items():
            # Find corresponding Kraken2 file
            kraken_fp = None
            if sample_id in kraken_files and mag_id in kraken_files[sample_id]:
                kraken_fp = kraken_files[sample_id][mag_id]
            
            tasks.append((mag_id, sample_id, mag_fp, kraken_fp, coverage_data, taxonomy_map))
    
    task_prep_time = time.time() - task_prep_start
    logger.info(f"Task preparation completed in {task_prep_time:.3f} seconds: {len(tasks)} MAGs to process")
    
    # Process all MAGs in parallel
    logger.info("Starting multiprocessing pool execution")
    parallel_start = time.time()
    results = {}
    processed_count = 0
    
    with multiprocessing.Pool(processes=n_cpus) as pool:
        for mag_id, sample_id, contigs_by_taxonomy in pool.map(_process_single_mag, tasks):
            processed_count += 1
            if sample_id not in results:
                results[sample_id] = []
            
            # Get BUSCO data for this MAG
            if mag_id in busco_results.index:
                logger.debug(f"Processing BUSCO data for MAG {mag_id}")
                busco_start = time.time()
                
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
                
                busco_time = time.time() - busco_start
                logger.debug(f"BUSCO data processing for MAG {mag_id} completed in {busco_time:.3f} seconds")
            else:
                logger.warning(f"No BUSCO data found for MAG {mag_id}, skipping")
    
    parallel_time = time.time() - parallel_start
    total_time = time.time() - start_time
    
    logger.info(f"Parallel processing completed in {parallel_time:.3f} seconds")
    logger.info(f"Total MAG processing completed in {total_time:.3f} seconds: {processed_count}/{len(tasks)} MAGs processed successfully")
    
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
    """Simplified main function with streamlined data processing."""
    start_time = time.time()
    logger.info("=== Starting MAG Quality Visualization ===")
    logger.info(f"Input data: {len(busco_results)} BUSCO results, {len(mag_taxonomy)} taxonomy entries, coverage matrix: {contig_coverage.shape}")
    
    # Prepare data
    logger.info("Preparing input data")
    prep_start = time.time()
    
    logger.debug("Converting mag_taxonomy to dictionary")
    taxonomy_map = mag_taxonomy.to_dict()
    
    logger.debug("Parsing coverage data")
    coverage_data = parse_coverage_data(contig_coverage)
    
    prep_time = time.time() - prep_start
    logger.info(f"Data preparation completed in {prep_time:.3f} seconds")
    
    # Process all MAGs in parallel
    logger.info("Starting MAG processing")
    processing_start = time.time()
    results = process_mags_parallel(
        mags, busco_results, kraken2_outputs, 
        taxonomy_map, coverage_data, n_cpus
    )
    processing_time = time.time() - processing_start
    logger.info(f"MAG processing completed in {processing_time:.3f} seconds")
    
    # Clean and save results
    logger.info("Cleaning and saving results")
    output_start = time.time()
    
    logger.debug("Cleaning data for JSON serialization")
    json_clean_start = time.time()
    cleaned_results = clean_for_json(results)
    json_clean_time = time.time() - json_clean_start
    logger.debug(f"JSON cleaning completed in {json_clean_time:.3f} seconds")
    
    # Count total results
    total_mags = sum(len(sample_results) for sample_results in cleaned_results.values())
    logger.info(f"Writing results for {len(cleaned_results)} samples with {total_mags} total MAGs")
    
    json_write_start = time.time()
    with open(os.path.join(output_dir, "data.json"), 'w') as f:
        json.dump(cleaned_results, f, indent=2)
    json_write_time = time.time() - json_write_start
    logger.debug(f"JSON file writing completed in {json_write_time:.3f} seconds")
    
    # Copy template files
    logger.debug("Copying template files")
    template_start = time.time()
    templates = [TEMPLATES / "mag_qc" / "index.html"]
    for fn in ("scripts.js", "bootstrapMagic.js", "styles.css"):
        shutil.copy(str(TEMPLATES / "mag_qc" / fn), os.path.join(output_dir, fn))
    
    q2templates.render(templates, output_dir, context={})
    template_time = time.time() - template_start
    logger.debug(f"Template processing completed in {template_time:.3f} seconds")
    
    output_time = time.time() - output_start
    total_time = time.time() - start_time
    
    logger.info(f"Output generation completed in {output_time:.3f} seconds")
    logger.info(f"=== MAG Quality Visualization completed in {total_time:.3f} seconds ===")
    
    # Summary logging
    logger.info("=== PERFORMANCE SUMMARY ===")
    logger.info(f"Data preparation: {prep_time:.3f}s ({prep_time/total_time*100:.1f}%)")
    logger.info(f"MAG processing: {processing_time:.3f}s ({processing_time/total_time*100:.1f}%)")
    logger.info(f"Output generation: {output_time:.3f}s ({output_time/total_time*100:.1f}%)")
    logger.info(f"Total runtime: {total_time:.3f}s")

