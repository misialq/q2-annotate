#!/usr/bin/env python3
"""
Script to parse MAG quality metrics and taxonomic assignments from BUSCO and Kraken2 results.

This script:
1. Parses BUSCO metrics file to get MAG quality data
2. Parses Kraken2 taxonomic assignment files organized by sample/MAG
3. Maps numeric taxon IDs to full taxonomy strings
4. Calculates Least Common Ancestor (LCA) for MAG taxonomy assignment
5. Outputs a comprehensive JSON file with all data
"""

import argparse
import json
import os
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from itertools import takewhile
from typing import List
import sys
import re


# Taxonomic ranks
RANKS = 'dkpcofgs'


def _find_lca(taxa):
    """Find least common ancestor between two semicolon-delimited strings.

    This code was copied over from https://github.com/bokulich-lab/RESCRIPt.
    """
    # LCA ends where zipped taxonomy strings no longer converge to len == 1
    taxa_comparison = [set(rank) - {None} for rank in zip(*taxa)]
    return (rank.pop() for rank in takewhile(
        lambda x: len(x) == 1, taxa_comparison))


def _fill_unclassified(row):
    if row.isnull().all():
        row[0] = 'Unclassified'
    return row


def _taxon_to_list(taxon_string, rank_handle=None):
    """Convert semicolon-delimited taxonomy string to list."""
    if pd.isna(taxon_string):
        return []
    
    # Split by semicolon
    taxa = taxon_string.split(';')
    
    # Remove rank handles if specified
    if rank_handle:
        taxa = [re.sub(rank_handle, '', taxon).strip() for taxon in taxa]
    
    # Remove empty strings and clean up
    taxa = [taxon.strip() for taxon in taxa if taxon.strip()]
    
    return taxa


def _join_ranks(taxa_list, ranks):
    """Join taxonomy list back to semicolon-delimited string with rank prefixes."""
    if not taxa_list:
        return 'Unclassified'
    
    result = []
    for i, taxon in enumerate(taxa_list):
        if i < len(ranks) and taxon:
            if not taxon.startswith(ranks[i]):
                result.append(f"{ranks[i]}{taxon}")
            else:
                result.append(taxon)
        elif taxon:
            result.append(taxon)
    
    return ';'.join(result) if result else 'Unclassified'


def _find_lcas(taxa_list: List[pd.DataFrame], mode: str = 'lca'):
    """Find the least common ancestor in every DataFrame of taxa.

    Args:
        taxa_list (List[pd.DataFrame]): A list of taxonomy DataFrames.
        mode (str): The mode used to determine the least common ancestor.

    Returns:
        pd.DataFrame: A DataFrame containing the LCA of each feature (MAG).
    """
    methods = {
        'lca': _find_lca,
        # 'super': _find_super_lca,
        # 'majority': _find_lca_majority
    }
    func = methods[mode]
    taxa = pd.concat(taxa_list)

    # Convert taxonomies to list; optionally remove rank handle
    taxa['Taxon'] = taxa['Taxon'].apply(
        lambda x: _taxon_to_list(x, rank_handle=f'^[{RANKS[:-1]}]__|s1?__')
    )

    # Find LCA for every MAG
    results = {}
    for mag_id in taxa['mag_id'].unique():
        data = taxa[taxa['mag_id'] == mag_id]['Taxon']
        result = list(func(data))
        results[mag_id] = result

    results = pd.DataFrame.from_dict(results, orient='index')
    results = results.apply(_fill_unclassified, axis=1)
    results = results.apply(lambda x: x.tolist(), axis=1).to_frame()
    results.columns = ['Taxon']

    # Join ranks
    ranks = [*[f'{r}__' for r in RANKS], 'ssp__']
    results['Taxon'] = results['Taxon'].apply(
        lambda x: _join_ranks(x, ranks)
    )

    results.index.name = 'Feature ID'
    return results


def parse_busco_file(busco_file):
    """Parse BUSCO results file to extract MAG quality metrics."""
    print(f"Parsing BUSCO file: {busco_file}")
    
    # Read BUSCO file
    df = pd.read_csv(busco_file, sep='\t')
    
    # Extract required columns
    busco_data = {}
    for _, row in df.iterrows():
        mag_id = row['mag_id']
        busco_data[mag_id] = {
            'sample_id': row['sample_id'],
            'completeness': row['completeness'],
            'contamination': row['contamination']
        }
    
    print(f"Loaded {len(busco_data)} MAG records from BUSCO file")
    return busco_data


def parse_taxonomy_file(taxonomy_file):
    """Parse taxonomy mapping file to map taxon IDs to full taxonomy."""
    print(f"Parsing taxonomy file: {taxonomy_file}")
    
    # Read taxonomy file
    df = pd.read_csv(taxonomy_file, sep='\t')
    
    # Create mapping from taxon ID to taxonomy
    taxonomy_map = {}
    for _, row in df.iterrows():
        taxon_id = row['Feature ID']
        taxonomy = row['Taxon']
        taxonomy_map[taxon_id] = taxonomy
    
    print(f"Loaded {len(taxonomy_map)} taxonomy mappings")
    return taxonomy_map


def extract_assigned_taxonomy_lca(mag_id, contigs_data, taxonomy_map):
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


def parse_kraken2_output(output_file, taxonomy_map):
    """Parse a single Kraken2 output file to get taxonomic assignments."""
    contigs_data = defaultdict(list)
    
    if not os.path.exists(output_file):
        print(f"Warning: Kraken2 output file not found: {output_file}")
        return contigs_data
    
    try:
        # Read Kraken2 output file
        with open(output_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                classification = parts[0]  # C = classified, U = unclassified
                contig_id = parts[1]       # Contig/sequence ID
                taxon_id = int(parts[2])   # Taxon ID
                length = int(parts[3])     # Sequence length
                
                # Only process classified sequences
                if classification == 'C' and taxon_id in taxonomy_map:
                    taxonomy = taxonomy_map[taxon_id]
                    
                    # Group by taxonomy at reasonable level (up to genus/species)
                    contigs_data[taxonomy].append({
                        'specific_taxonomy': taxonomy,
                        'length': length
                    })
        
    except Exception as e:
        print(f"Error parsing {output_file}: {e}")
    
    return dict(contigs_data)


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


def determine_quality(completeness, contamination):
    """Determine MAG quality based on completeness and contamination."""
    if completeness >= 90 and contamination <= 5:
        return "High"
    elif completeness >= 50 and contamination <= 10:
        return "Medium"
    else:
        return "Low"


def parse_kraken2_directory(kraken2_dir, taxonomy_map):
    """Parse directory structure containing Kraken2 output files."""
    print(f"Parsing Kraken2 directory: {kraken2_dir}")
    
    mag_taxonomies = {}
    
    # Walk through directory structure: sample_id/mag_id.output.txt
    for sample_dir in os.listdir(kraken2_dir):
        sample_path = os.path.join(kraken2_dir, sample_dir)
        
        if not os.path.isdir(sample_path):
            continue
        
        print(f"Processing sample: {sample_dir}")
        
        for filename in os.listdir(sample_path):
            if filename.endswith('.output.txt'):
                # Extract MAG ID from filename (remove .output.txt extension)
                mag_id = filename.replace('.output.txt', '')
                output_file = os.path.join(sample_path, filename)
                
                # Parse Kraken2 output for this MAG
                contigs_data = parse_kraken2_output(output_file, taxonomy_map)
                
                # Store the results
                mag_taxonomies[mag_id] = {
                    'sample_id': sample_dir,
                    'contigs_data': contigs_data
                }
    
    print(f"Processed taxonomic data for {len(mag_taxonomies)} MAGs")
    return mag_taxonomies


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


def main():
    parser = argparse.ArgumentParser(description='Parse MAG quality and taxonomic data')
    parser.add_argument('--busco', required=True, help='Path to BUSCO results TSV file')
    parser.add_argument('--kraken2-dir', required=True, help='Path to directory containing Kraken2 output files (sample_id/mag_id.output.txt)')
    parser.add_argument('--taxonomy', required=True, help='Path to taxonomy mapping TSV file')
    parser.add_argument('--output', required=True, help='Path to output JSON file')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.busco):
        print(f"Error: BUSCO file not found: {args.busco}")
        sys.exit(1)
    
    if not os.path.exists(args.kraken2_dir):
        print(f"Error: Kraken2 directory not found: {args.kraken2_dir}")
        sys.exit(1)
    
    if not os.path.exists(args.taxonomy):
        print(f"Error: Taxonomy file not found: {args.taxonomy}")
        sys.exit(1)
    
    try:
        # Parse input files
        busco_data = parse_busco_file(args.busco)
        taxonomy_map = parse_taxonomy_file(args.taxonomy)
        mag_taxonomies = parse_kraken2_directory(args.kraken2_dir, taxonomy_map)
        
        # Combine data and organize by sample
        result = defaultdict(list)
        
        for mag_id, busco_info in busco_data.items():
            sample_id = busco_info['sample_id']
            
            # Get taxonomic data for this MAG
            taxonomy_info = mag_taxonomies.get(mag_id, {
                'sample_id': sample_id,
                'contigs_data': {}
            })
            
            contigs_data = taxonomy_info['contigs_data']
            
            # Calculate MAG statistics
            stats = calculate_mag_stats(contigs_data)
            
            # Determine assigned taxonomy using LCA
            assigned_taxonomy = extract_assigned_taxonomy_lca(mag_id, contigs_data, taxonomy_map)
            
            # Determine quality
            quality = determine_quality(busco_info['completeness'], busco_info['contamination'])
            
            # Create MAG entry
            mag_entry = {
                'id': mag_id,
                'contamination': busco_info['contamination'],
                'completeness': busco_info['completeness'],
                'assigned_taxonomy': assigned_taxonomy,
                'contigs': contigs_data,
                'n50': stats['n50'],
                'l50': stats['l50'],
                'total_length': stats['total_length'],
                'quality': quality
            }
            
            result[sample_id].append(mag_entry)
        
        # Write output JSON
        print(f"Writing output to: {args.output}")
        
        # Clean NaN values before JSON serialization
        cleaned_result = clean_nan_values(dict(result))
        
        with open(args.output, 'w') as f:
            json.dump(cleaned_result, f, indent=2)
        
        print(f"Successfully processed {len(busco_data)} MAGs across {len(result)} samples")
        
        # Print summary
        total_mags = sum(len(mags) for mags in result.values())
        print(f"Summary:")
        print(f"  Total samples: {len(result)}")
        print(f"  Total MAGs: {total_mags}")
        
        for sample_id, mags in result.items():
            print(f"  {sample_id}: {len(mags)} MAGs")
    
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main() 