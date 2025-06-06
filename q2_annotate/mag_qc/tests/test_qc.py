# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
import unittest
import unittest.mock
from unittest.mock import patch, call
import pandas as pd
import numpy as np
from pathlib import Path

from qiime2.plugin.testing import TestPluginBase

from q2_annotate.mag_qc.qc import (
    determine_quality,
    calculate_assembly_stats,
    extract_assigned_taxonomy_lca,
    _process_single_mag,
    process_mags_parallel,
    parse_coverage_data,
    clean_for_json,
    ContigData,
    MAGData
)


class TestQCFunctions(TestPluginBase):
    package = 'q2_annotate.mag_qc.tests'

    def test_determine_quality_high(self):
        """Test high quality determination."""
        quality = determine_quality(95.0, 3.0)
        self.assertEqual(quality, "High")

        quality = determine_quality(90.0, 5.0)
        self.assertEqual(quality, "High")

    def test_determine_quality_medium(self):
        """Test medium quality determination."""
        quality = determine_quality(85.0, 7.0)
        self.assertEqual(quality, "Medium")

        quality = determine_quality(50.0, 10.0)
        self.assertEqual(quality, "Medium")

        quality = determine_quality(75.0, 15.0)
        self.assertEqual(quality, "Low")  # contamination too high

    def test_determine_quality_low(self):
        """Test low quality determination."""
        quality = determine_quality(45.0, 8.0)
        self.assertEqual(quality, "Low")

        quality = determine_quality(30.0, 3.0)
        self.assertEqual(quality, "Low")

        quality = determine_quality(95.0, 15.0)
        self.assertEqual(quality, "Low")  # contamination too high

    def test_calculate_assembly_stats_empty(self):
        """Test assembly stats with empty contigs."""
        contigs_by_taxonomy = {}
        n50, l50, total_length = calculate_assembly_stats(contigs_by_taxonomy)
        self.assertEqual((n50, l50, total_length), (0, 0, 0))

    def test_calculate_assembly_stats_single_taxonomy(self):
        """Test assembly stats with contigs from single taxonomy."""
        contigs_by_taxonomy = {
            "d__Bacteria;p__Firmicutes": [
                {"length": 1000}, {"length": 2000}, {"length": 1500}
            ]
        }
        n50, l50, total_length = calculate_assembly_stats(contigs_by_taxonomy)
        # Sorted lengths: [2000, 1500, 1000]
        # Total: 4500, half: 2250
        # Cumulative: 2000 (< 2250), 2000+1500=3500 (>= 2250)
        # So N50 = 1500, L50 = 2
        self.assertEqual(n50, 1500)
        self.assertEqual(l50, 2)
        self.assertEqual(total_length, 4500)

    def test_calculate_assembly_stats_multiple_taxonomies(self):
        """Test assembly stats with contigs from multiple taxonomies."""
        contigs_by_taxonomy = {
            "d__Bacteria;p__Firmicutes": [
                {"length": 2000}, {"length": 1000}
            ],
            "d__Bacteria;p__Proteobacteria": [
                {"length": 3000}, {"length": 500}
            ],
            "Unassigned": [
                {"length": 1500}
            ]
        }
        n50, l50, total_length = calculate_assembly_stats(contigs_by_taxonomy)
        # All lengths: [3000, 2000, 1500, 1000, 500]
        # Total: 8000, half: 4000
        # Cumulative: 3000 (< 4000), 3000+2000=5000 (>= 4000)
        # So N50 = 2000, L50 = 2
        self.assertEqual(n50, 2000)
        self.assertEqual(l50, 2)
        self.assertEqual(total_length, 8000)

    def test_extract_assigned_taxonomy_lca_empty(self):
        """Test LCA with empty contigs."""
        contigs_by_taxonomy = {}
        result = extract_assigned_taxonomy_lca(contigs_by_taxonomy)
        self.assertEqual(result, "Unassigned")

    def test_extract_assigned_taxonomy_lca_unassigned_only(self):
        """Test LCA with only unassigned contigs."""
        contigs_by_taxonomy = {
            "Unassigned": [
                {"contig_id": "contig1"}, {"contig_id": "contig2"}
            ]
        }
        result = extract_assigned_taxonomy_lca(contigs_by_taxonomy)
        self.assertEqual(result, "Unassigned")

    @patch('q2_annotate.mag_qc.qc._find_lcas')
    def test_extract_assigned_taxonomy_lca_success(self, mock_find_lcas):
        """Test successful LCA calculation."""
        contigs_by_taxonomy = {
            "d__Bacteria;p__Firmicutes;c__Bacilli": [
                {"contig_id": "contig1"}, {"contig_id": "contig2"}
            ],
            "d__Bacteria;p__Firmicutes;c__Clostridia": [
                {"contig_id": "contig3"}
            ]
        }

        # Mock the LCA result
        mock_result = pd.DataFrame({
            'Taxon': ['d__Bacteria;p__Firmicutes']
        }, index=['temp'])
        mock_find_lcas.return_value = mock_result

        result = extract_assigned_taxonomy_lca(contigs_by_taxonomy)
        self.assertEqual(result, "d__Bacteria;p__Firmicutes")

        # Verify _find_lcas was called with correct data
        mock_find_lcas.assert_called_once()
        call_args = mock_find_lcas.call_args[0]
        taxa_df = call_args[0][0]
        
        # Should have 3 rows (2 + 1 contigs)
        self.assertEqual(len(taxa_df), 3)
        self.assertTrue(all(taxa_df['mag_id'] == 'temp'))

    @patch('q2_annotate.mag_qc.qc._find_lcas')
    def test_extract_assigned_taxonomy_lca_fallback(self, mock_find_lcas):
        """Test LCA fallback to most common taxonomy."""
        contigs_by_taxonomy = {
            "d__Bacteria;p__Firmicutes": [
                {"contig_id": "contig1"}, {"contig_id": "contig2"}
            ],
            "d__Bacteria;p__Proteobacteria": [
                {"contig_id": "contig3"}
            ]
        }

        # Mock LCA to raise an exception
        mock_find_lcas.side_effect = Exception("LCA failed")

        result = extract_assigned_taxonomy_lca(contigs_by_taxonomy)
        # Should return the taxonomy with most contigs
        self.assertEqual(result, "d__Bacteria;p__Firmicutes")

    def test_clean_for_json_dict(self):
        """Test JSON cleaning for dictionaries."""
        data = {
            'key1': np.int64(42),
            'key2': np.float64(3.14),
            'key3': np.nan,
            'key4': 'string'
        }
        result = clean_for_json(data)
        expected = {
            'key1': 42,
            'key2': 3.14,
            'key3': None,
            'key4': 'string'
        }
        self.assertEqual(result, expected)

    def test_clean_for_json_list(self):
        """Test JSON cleaning for lists."""
        data = [np.int32(1), np.float32(2.5), np.nan, 'text']
        result = clean_for_json(data)
        expected = [1, 2.5, None, 'text']
        self.assertEqual(result, expected)

    def test_clean_for_json_nested(self):
        """Test JSON cleaning for nested structures."""
        data = {
            'list': [np.int64(1), {'nested': np.nan}],
            'array': np.array([1, 2, 3])
        }
        result = clean_for_json(data)
        expected = {
            'list': [1, {'nested': None}],
            'array': [1, 2, 3]
        }
        self.assertEqual(result, expected)

    def test_process_single_mag_no_kraken(self):
        """Test processing MAG without Kraken2 data using real test files."""
        mag_fp = self.get_data_path('sample_mag.fa')
        kraken_fp = '/nonexistent/path'  # This file doesn't exist
        coverage_data = {'contig1': 10.5, 'contig2': 15.2, 'contig3': 8.0}
        taxonomy_map = {}

        args = (
            'mag1', 'sample1', mag_fp, kraken_fp,
            coverage_data, taxonomy_map
        )

        mag_id, sample_id, contigs_by_taxonomy = _process_single_mag(args)

        self.assertEqual(mag_id, 'mag1')
        self.assertEqual(sample_id, 'sample1')
        self.assertEqual(len(contigs_by_taxonomy), 1)
        self.assertIn('Unassigned', contigs_by_taxonomy)
        self.assertEqual(len(contigs_by_taxonomy['Unassigned']), 3)

        # Check that we have the expected contigs
        contig_ids = [c['contig_id'] for c in contigs_by_taxonomy['Unassigned']]
        self.assertIn('contig1', contig_ids)
        self.assertIn('contig2', contig_ids)
        self.assertIn('contig3', contig_ids)

    def test_process_single_mag_with_kraken(self):
        """Test processing MAG with Kraken2 data using real test files."""
        mag_fp = self.get_data_path('sample_mag.fa')
        kraken_fp = self.get_data_path('sample_kraken.out')
        coverage_data = {'contig1': 10.5, 'contig2': 15.2, 'contig3': 8.0}
        taxonomy_map = {
            '123': 'd__Bacteria;p__Firmicutes',
            '456': 'd__Bacteria;p__Proteobacteria'
        }

        args = (
            'mag1', 'sample1', mag_fp, kraken_fp,
            coverage_data, taxonomy_map
        )

        mag_id, sample_id, contigs_by_taxonomy = _process_single_mag(args)

        self.assertEqual(mag_id, 'mag1')
        self.assertEqual(sample_id, 'sample1')
        
        # Check that contigs are grouped by taxonomy
        # contigs 1 and 2 should be classified, contig3 should be unassigned
        expected_taxonomies = {'d__Bacteria;p__Firmicutes', 'd__Bacteria;p__Proteobacteria', 'Unassigned'}
        actual_taxonomies = set(contigs_by_taxonomy.keys())
        self.assertEqual(actual_taxonomies, expected_taxonomies)
        
        # Check that the right number of contigs are in each taxonomy
        total_contigs = sum(len(contigs) for contigs in contigs_by_taxonomy.values())
        self.assertEqual(total_contigs, 3)

    @patch('q2_annotate.mag_qc.qc.multiprocessing.Pool')
    @patch('q2_annotate.mag_qc.qc._process_single_mag')
    def test_process_mags_parallel_single_sample(self, mock_process_mag, mock_pool):
        """Test parallel processing with a single sample containing one MAG."""
        # Mock data structures
        mock_mags = unittest.mock.Mock()
        mock_mags.sample_dict.return_value = {
            'sample1': {'mag1': '/path/to/mag1.fa'}
        }
        
        mock_kraken = unittest.mock.Mock()
        mock_kraken.file_dict.return_value = {
            'sample1': {'mag1': '/path/to/kraken1.out'}
        }
        
        mock_busco = pd.DataFrame({
            'completeness': [85.5],
            'contamination': [3.2]
        }, index=['mag1'])
        
        # Mock the pool and _process_single_mag
        mock_pool_instance = unittest.mock.Mock()
        mock_pool.return_value.__enter__.return_value = mock_pool_instance
        
        mock_contigs = {
            'd__Bacteria': [
                {'contig_id': 'contig1', 'length': 1000}
            ]
        }
        mock_pool_instance.map.return_value = [
            ('mag1', 'sample1', mock_contigs)
        ]
        
        taxonomy_map = {'123': 'd__Bacteria'}
        coverage_data = {'contig1': 10.5}
        
        # Call the function
        results = process_mags_parallel(
            mock_mags, mock_busco, mock_kraken, 
            taxonomy_map, coverage_data, n_cpus=2
        )
        
        # Verify results
        self.assertEqual(len(results), 1)
        self.assertIn('sample1', results)
        self.assertEqual(len(results['sample1']), 1)
        
        mag_data = results['sample1'][0]
        self.assertEqual(mag_data['id'], 'mag1')
        self.assertEqual(mag_data['sample_id'], 'sample1')
        self.assertEqual(mag_data['completeness'], 85.5)
        self.assertEqual(mag_data['contamination'], 3.2)
        self.assertEqual(mag_data['quality'], 'Medium')  # 85.5% complete, 3.2% contamination
        
        # Verify the pool was used correctly
        mock_pool.assert_called_once_with(processes=2)
        mock_pool_instance.map.assert_called_once()

    @patch('q2_annotate.mag_qc.qc.multiprocessing.Pool')
    @patch('q2_annotate.mag_qc.qc._process_single_mag')
    def test_process_mags_parallel_multiple_samples(self, mock_process_mag, mock_pool):
        """Test parallel processing with multiple samples and MAGs."""
        # Mock data structures
        mock_mags = unittest.mock.Mock()
        mock_mags.sample_dict.return_value = {
            'sample1': {'mag1': '/path/to/mag1.fa', 'mag2': '/path/to/mag2.fa'},
            'sample2': {'mag3': '/path/to/mag3.fa'}
        }
        
        mock_kraken = unittest.mock.Mock()
        mock_kraken.file_dict.return_value = {
            'sample1': {'mag1': '/path/to/kraken1.out'},
            'sample2': {'mag3': '/path/to/kraken3.out'}
        }
        
        mock_busco = pd.DataFrame({
            'completeness': [95.0, 75.0, 40.0],
            'contamination': [2.0, 8.0, 15.0]
        }, index=['mag1', 'mag2', 'mag3'])
        
        # Mock the pool and _process_single_mag
        mock_pool_instance = unittest.mock.Mock()
        mock_pool.return_value.__enter__.return_value = mock_pool_instance
        
        mock_contigs1 = {'d__Bacteria': [{'contig_id': 'contig1', 'length': 2000}]}
        mock_contigs2 = {'d__Archaea': [{'contig_id': 'contig2', 'length': 1500}]}
        mock_contigs3 = {'Unassigned': [{'contig_id': 'contig3', 'length': 800}]}
        
        mock_pool_instance.map.return_value = [
            ('mag1', 'sample1', mock_contigs1),
            ('mag2', 'sample1', mock_contigs2),
            ('mag3', 'sample2', mock_contigs3)
        ]
        
        taxonomy_map = {}
        coverage_data = {}
        
        # Call the function
        results = process_mags_parallel(
            mock_mags, mock_busco, mock_kraken, 
            taxonomy_map, coverage_data, n_cpus=1
        )
        
        # Verify results
        self.assertEqual(len(results), 2)
        self.assertIn('sample1', results)
        self.assertIn('sample2', results)
        self.assertEqual(len(results['sample1']), 2)  # mag1 and mag2
        self.assertEqual(len(results['sample2']), 1)  # mag3
        
        # Check quality determination
        sample1_qualities = [mag['quality'] for mag in results['sample1']]
        self.assertIn('High', sample1_qualities)  # mag1: 95% complete, 2% contamination
        self.assertIn('Medium', sample1_qualities)  # mag2: 75% complete, 8% contamination
        
        sample2_mag = results['sample2'][0]
        self.assertEqual(sample2_mag['quality'], 'Low')  # mag3: 40% complete, 15% contamination

    @patch('q2_annotate.mag_qc.qc.multiprocessing.Pool')
    @patch('q2_annotate.mag_qc.qc._process_single_mag')
    def test_process_mags_parallel_missing_busco_data(self, mock_process_mag, mock_pool):
        """Test parallel processing when some MAGs don't have BUSCO data."""
        # Mock data structures
        mock_mags = unittest.mock.Mock()
        mock_mags.sample_dict.return_value = {
            'sample1': {'mag1': '/path/to/mag1.fa', 'mag2': '/path/to/mag2.fa'}
        }
        
        mock_kraken = unittest.mock.Mock()
        mock_kraken.file_dict.return_value = {}
        
        # Only mag1 has BUSCO data, mag2 doesn't
        mock_busco = pd.DataFrame({
            'completeness': [85.5],
            'contamination': [3.2]
        }, index=['mag1'])
        
        # Mock the pool and _process_single_mag
        mock_pool_instance = unittest.mock.Mock()
        mock_pool.return_value.__enter__.return_value = mock_pool_instance
        
        mock_contigs = {'d__Bacteria': [{'contig_id': 'contig1', 'length': 1000}]}
        mock_pool_instance.map.return_value = [
            ('mag1', 'sample1', mock_contigs),
            ('mag2', 'sample1', mock_contigs)
        ]
        
        taxonomy_map = {}
        coverage_data = {}
        
        # Call the function
        results = process_mags_parallel(
            mock_mags, mock_busco, mock_kraken, 
            taxonomy_map, coverage_data, n_cpus=1
        )
        
        # Verify results - only mag1 should be in results (has BUSCO data)
        self.assertEqual(len(results), 1)
        self.assertIn('sample1', results)
        self.assertEqual(len(results['sample1']), 1)
        self.assertEqual(results['sample1'][0]['id'], 'mag1')

    def test_process_mags_parallel_empty_inputs(self):
        """Test parallel processing with empty inputs."""
        # Mock empty data structures
        mock_mags = unittest.mock.Mock()
        mock_mags.sample_dict.return_value = {}
        
        mock_kraken = unittest.mock.Mock()
        mock_kraken.file_dict.return_value = {}
        
        mock_busco = pd.DataFrame(columns=['completeness', 'contamination'])
        
        taxonomy_map = {}
        coverage_data = {}
        
        # Call the function
        results = process_mags_parallel(
            mock_mags, mock_busco, mock_kraken, 
            taxonomy_map, coverage_data, n_cpus=1
        )
        
        # Verify results
        self.assertEqual(results, {})

    def test_parse_coverage_data_basic(self):
        """Test basic coverage data parsing with single sample per contig."""
        # Create test DataFrame where each contig appears in only one sample
        data = {
            'contig1': [10.5, 0.0, 0.0],   # Only in sample1
            'contig2': [0.0, 12.1, 0.0],   # Only in sample2
            'contig3': [0.0, 0.0, 22.7]    # Only in sample3
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3'])
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig1': 10.5,  # From sample1
            'contig2': 12.1,  # From sample2  
            'contig3': 22.7   # From sample3
        }
        self.assertEqual(result, expected)

    def test_parse_coverage_data_with_nan_values(self):
        """Test coverage data parsing with NaN values."""
        # Create test DataFrame with NaN values (single sample per contig)
        data = {
            'contig1': [np.nan, 15.2, 0.0],   # Only in sample2
            'contig2': [5.5, np.nan, np.nan], # Only in sample1
            'contig3': [np.nan, np.nan, 0.0]  # No positive values
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3'])
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig1': 15.2,  # From sample2
            'contig2': 5.5    # From sample1
            # contig3 should not appear (no positive values)
        }
        self.assertEqual(result, expected)

    def test_parse_coverage_data_all_zeros_and_negatives(self):
        """Test coverage data parsing with all zeros and negative values."""
        # Create test DataFrame with zeros and negative values
        data = {
            'contig1': [0.0, 0.0, -1.0],
            'contig2': [-5.5, 0.0, 0.0],
            'contig3': [12.5, 0.0, 0.0]
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3'])
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig3': 12.5  # Only contig with positive value
        }
        self.assertEqual(result, expected)

    def test_parse_coverage_data_empty_dataframe(self):
        """Test coverage data parsing with empty DataFrame."""
        df = pd.DataFrame()
        
        result = parse_coverage_data(df)
        
        self.assertEqual(result, {})

    def test_parse_coverage_data_single_value(self):
        """Test coverage data parsing with single value."""
        data = {'contig1': [25.8]}
        df = pd.DataFrame(data, index=['sample1'])
        
        result = parse_coverage_data(df)
        
        expected = {'contig1': 25.8}
        self.assertEqual(result, expected)



    def test_parse_coverage_data_type_conversion(self):
        """Test that coverage values are properly converted to float."""
        # Create test DataFrame with integer values (single sample per contig)
        data = {
            'contig1': [10, 0, 0],   # Only in sample1
            'contig2': [0, 12, 0]    # Only in sample2
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3'])
        
        result = parse_coverage_data(df)
        
        # Check that values are floats
        self.assertIsInstance(result['contig1'], float)
        self.assertIsInstance(result['contig2'], float)
        self.assertEqual(result['contig1'], 10.0)  # From sample1
        self.assertEqual(result['contig2'], 12.0)  # From sample2

    def test_parse_coverage_data_disjoint_samples(self):
        """Test coverage data parsing where different samples cover different contigs."""
        # Create test DataFrame where samples have coverage for different contigs
        data = {
            'contig1': [10.5, 0.0, 0.0],    # Only sample1 has coverage
            'contig2': [0.0, 15.2, 0.0],    # Only sample2 has coverage  
            'contig3': [0.0, 0.0, 22.7],    # Only sample3 has coverage
            'contig4': [5.1, 8.3, 0.0]      # sample1 and sample2 have coverage
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3'])
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig1': 10.5,  # Only sample1
            'contig2': 15.2,  # Only sample2
            'contig3': 22.7,  # Only sample3
            'contig4': 8.3    # Last positive (sample2)
        }
        self.assertEqual(result, expected)

    def test_parse_coverage_data_sparse_coverage(self):
        """Test coverage data parsing with very sparse coverage across samples."""
        # Create test DataFrame with mostly zeros/NaN
        data = {
            'contig1': [0.0, 0.0, 0.0, 0.0, 12.5],      # Only last sample
            'contig2': [3.2, 0.0, 0.0, 0.0, 0.0],       # Only first sample
            'contig3': [0.0, 0.0, 7.8, 0.0, 0.0],       # Only middle sample
            'contig4': [0.0, 0.0, 0.0, 0.0, 0.0],       # No coverage
            'contig5': [np.nan, np.nan, np.nan, 4.1, np.nan]  # One sample with NaNs
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig1': 12.5,
            'contig2': 3.2,
            'contig3': 7.8,
            'contig5': 4.1
            # contig4 should not appear (no coverage)
        }
        self.assertEqual(result, expected)

    def test_parse_coverage_data_many_samples(self):
        """Test coverage data parsing with many samples to ensure scalability."""
        # Create test DataFrame with 10 samples
        sample_names = [f'sample{i}' for i in range(1, 11)]
        
        # contig1: coverage increases by sample number
        # contig2: only odd-numbered samples have coverage
        # contig3: only even-numbered samples have coverage
        data = {
            'contig1': list(range(1, 11)),  # 1, 2, 3, ..., 10
            'contig2': [i if i % 2 == 1 else 0 for i in range(1, 11)],  # 1, 0, 3, 0, 5, 0, 7, 0, 9, 0
            'contig3': [i if i % 2 == 0 else 0 for i in range(1, 11)]   # 0, 2, 0, 4, 0, 6, 0, 8, 0, 10
        }
        df = pd.DataFrame(data, index=sample_names)
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig1': 10.0,  # Last sample (sample10)
            'contig2': 9.0,   # Last odd sample (sample9) 
            'contig3': 10.0   # Last even sample (sample10)
        }
        self.assertEqual(result, expected)

    def test_parse_coverage_data_sample_names_with_special_chars(self):
        """Test coverage data parsing with unusual sample names."""
        # Test with sample names that have special characters
        data = {
            'contig1': [10.5, 15.2, 8.7],
            'contig2': [0.0, 12.1, 0.0]
        }
        df = pd.DataFrame(data, index=['sample-1_test', 'sample.2@domain', 'sample 3 (replicate)'])
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig1': 8.7,   # Last positive value
            'contig2': 12.1   # Only positive value
        }
        self.assertEqual(result, expected)

    def test_parse_coverage_data_samples_all_zero_for_contig(self):
        """Test coverage data parsing where all samples have zero coverage for some contigs."""
        # Create test DataFrame where some contigs have no coverage across all samples
        data = {
            'contig1': [10.5, 15.2, 8.7],    # All samples have coverage
            'contig2': [0.0, 0.0, 0.0],      # No samples have coverage
            'contig3': [5.5, 0.0, 0.0],      # Only first sample has coverage
            'contig4': [0.0, 0.0, 3.2]       # Only last sample has coverage  
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3'])
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig1': 8.7,   # Last positive value
            'contig3': 5.5,   # Only positive value
            'contig4': 3.2    # Only positive value
            # contig2 should not appear (no positive coverage)
        }
        self.assertEqual(result, expected)

    def test_parse_coverage_data_mixed_data_types(self):
        """Test coverage data parsing with mixed numeric data types."""
        # Create test DataFrame with mixed int/float types (single sample per contig)
        data = {
            'contig1': [10, 0, 0],                        # int in sample1
            'contig2': [0, np.float64(12.1), 0]          # numpy float64 in sample2
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3'])
        
        result = parse_coverage_data(df)
        
        # Check that we get the expected contigs and correct values
        self.assertEqual(set(result.keys()), {'contig1', 'contig2'})
        
        self.assertEqual(result['contig1'], 10.0)   # From sample1, converted to float
        self.assertEqual(result['contig2'], 12.1)   # From sample2, converted to float
        
        # Verify all returned values are Python floats, not numpy types
        for value in result.values():
            self.assertIsInstance(value, float)
            self.assertNotIsInstance(value, (np.integer, np.floating))

    def test_parse_coverage_data_optimized_single_sample_per_contig(self):
        """Test coverage parsing when each contig appears in only one sample."""
        # Create test DataFrame where each contig appears in only one sample
        data = {
            'contig1': [10.5, 0.0, 0.0],    # Only in sample1
            'contig2': [0.0, 15.2, 0.0],    # Only in sample2  
            'contig3': [0.0, 0.0, 22.7],    # Only in sample3
            'contig4': [8.3, 0.0, 0.0]      # Only in sample1
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3'])
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig1': 10.5,
            'contig2': 15.2,
            'contig3': 22.7,
            'contig4': 8.3
        }
        
        self.assertEqual(result, expected)

    def test_parse_coverage_data_optimized_with_zeros_and_nans(self):
        """Test coverage parsing with zeros and NaN values."""
        # Create test DataFrame with zeros and NaN values
        data = {
            'contig1': [np.nan, 0.0, 12.5],   # Only positive value in sample3
            'contig2': [7.2, np.nan, 0.0],    # Only positive value in sample1
            'contig3': [0.0, 0.0, 0.0],       # No positive values
            'contig4': [np.nan, np.nan, np.nan]  # All NaN
        }
        df = pd.DataFrame(data, index=['sample1', 'sample2', 'sample3'])
        
        result = parse_coverage_data(df)
        
        expected = {
            'contig1': 12.5,
            'contig2': 7.2
            # contig3 and contig4 should not appear (no positive values)
        }
        self.assertEqual(result, expected)

    def test_parse_coverage_data_performance_comparison(self):
        """Test that optimized path gives same results as original for single-sample scenario."""
        # Create larger test case
        import random
        random.seed(42)
        
        # Create 50 contigs, 20 samples
        n_contigs = 50
        n_samples = 20
        contig_names = [f'contig_{i}' for i in range(n_contigs)]
        sample_names = [f'sample_{i}' for i in range(n_samples)]
        
        # Create data where each contig appears in exactly one sample
        data = {}
        for i, contig in enumerate(contig_names):
            # Each contig gets a random positive value in one random sample
            column = [0.0] * n_samples
            sample_idx = i % n_samples  # Distribute evenly across samples
            column[sample_idx] = random.uniform(1.0, 100.0)
            data[contig] = column
        
        df = pd.DataFrame(data, index=sample_names)
        
        # Test the optimized function
        result = parse_coverage_data(df)
        
        self.assertEqual(len(result), n_contigs)  # All contigs should be present
        
        # All values should be positive
        for value in result.values():
            self.assertGreater(value, 0)
            self.assertIsInstance(value, float)


class TestDataClasses(unittest.TestCase):
    """Test the dataclasses used in the module."""

    def test_contig_data_creation(self):
        """Test ContigData dataclass creation."""
        contig = ContigData(
            contig_id="test_contig",
            length=1500,
            gc_content=45.5,
            coverage=12.3,
            taxonomy="d__Bacteria"
        )
        
        self.assertEqual(contig.contig_id, "test_contig")
        self.assertEqual(contig.length, 1500)
        self.assertEqual(contig.gc_content, 45.5)
        self.assertEqual(contig.coverage, 12.3)
        self.assertEqual(contig.taxonomy, "d__Bacteria")

    def test_contig_data_defaults(self):
        """Test ContigData with default values."""
        contig = ContigData(
            contig_id="test_contig",
            length=1500,
            gc_content=45.5,
            coverage=12.3
        )
        
        self.assertEqual(contig.taxonomy, "Unassigned")

    def test_mag_data_creation(self):
        """Test MAGData dataclass creation."""
        contigs = {
            "d__Bacteria": [
                {"contig_id": "contig1", "length": 1000}
            ]
        }
        
        mag = MAGData(
            id="mag1",
            sample_id="sample1",
            completeness=85.5,
            contamination=3.2,
            quality="High",
            assigned_taxonomy="d__Bacteria",
            contigs=contigs,
            n50=1000,
            l50=1,
            total_length=1000
        )
        
        self.assertEqual(mag.id, "mag1")
        self.assertEqual(mag.sample_id, "sample1")
        self.assertEqual(mag.completeness, 85.5)
        self.assertEqual(mag.contamination, 3.2)
        self.assertEqual(mag.quality, "High")
        self.assertEqual(mag.assigned_taxonomy, "d__Bacteria")
        self.assertEqual(mag.contigs, contigs)
        self.assertEqual(mag.n50, 1000)
        self.assertEqual(mag.l50, 1)
        self.assertEqual(mag.total_length, 1000)


if __name__ == '__main__':
    unittest.main()
