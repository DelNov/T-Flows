import numpy as np
import sys

def read_cfd_file(filename):
    """
    Read CFD output file with format:
    cell_id value
    """
    try:
        data = np.loadtxt(filename)
        # If file has only one column (just values without cell IDs)
        if data.ndim == 1:
            cell_ids = np.arange(1, len(data) + 1)
            values = data
        else:
            cell_ids = data[:, 0].astype(int)
            values = data[:, 1]
        return cell_ids, values
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        return None, None

def compare_cfd_files(file1, file2):
    """
    Compare two CFD files and return statistics
    """
    print(f"Comparing files:")
    print(f"  Reference: {file1}")
    print(f"  Test: {file2}")
    print()
    
    # Read files
    cell_ids1, values1 = read_cfd_file(file1)
    cell_ids2, values2 = read_cfd_file(file2)
    
    if values1 is None or values2 is None:
        print("Error: Could not read one or both files")
        return
    
    # Check if files have same number of data points
    if len(values1) != len(values2):
        print(f"Warning: Files have different lengths ({len(values1)} vs {len(values2)})")
        min_len = min(len(values1), len(values2))
        values1 = values1[:min_len]
        values2 = values2[:min_len]
        cell_ids1 = cell_ids1[:min_len]
        print(f"Comparing first {min_len} cells only")
    
    # Calculate differences
    differences = np.abs(values1 - values2)
    
    # Check if files are identical
    max_diff = np.max(differences)
    if max_diff == 0:
        print("FILES ARE IDENTICAL")
        return 'identical', None, None, None
    
    # Avoid division by zero for relative differences
    with np.errstate(divide='ignore', invalid='ignore'):
        relative_diff = np.abs(values1 - values2) / (np.abs(values1))
        relative_diff = np.nan_to_num(relative_diff, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Statistics
    max_diff_idx = np.argmax(differences)
    max_relative_idx = np.argmax(relative_diff)
    
    stats = {
        'mean_absolute_diff': np.mean(differences),
        'max_absolute_diff': np.max(differences),
        'mean_relative_diff': np.mean(relative_diff),
        'max_relative_diff': np.max(relative_diff),
        'cell_max_abs_diff': cell_ids1[max_diff_idx],
        'cell_max_rel_diff': cell_ids1[max_relative_idx],
        'value1_at_max_abs': values1[max_diff_idx],
        'value2_at_max_abs': values2[max_diff_idx],
        'value1_at_max_rel': values1[max_relative_idx],
        'value2_at_max_rel': values2[max_relative_idx],
        'rms_diff': np.sqrt(np.mean(differences**2)),
        'total_cells': len(values1)
    }
    
    return stats, differences, relative_diff, cell_ids1, values1, values2

def print_comparison_results(stats):
    """Print formatted comparison results"""
    print("=" * 60)
    print("CFD FILE COMPARISON RESULTS")
    print("=" * 60)
    print(f"Total cells compared: {stats['total_cells']}")
    print()
    print("ABSOLUTE DIFFERENCES:")
    print(f"  Mean absolute difference: {stats['mean_absolute_diff']:.6e}")
    print(f"  Maximum absolute difference: {stats['max_absolute_diff']:.6e}")
    print(f"  RMS difference: {stats['rms_diff']:.6e}")
    print()
    print("RELATIVE DIFFERENCES:")
    print(f"  Mean relative difference: {stats['mean_relative_diff']:.6e}")
    print(f"  Maximum relative difference: {stats['max_relative_diff']:.6e}")
    print()
    print("LOCATIONS OF MAXIMUM DIFFERENCES:")
    print(f"  Cell with max absolute difference: {stats['cell_max_abs_diff']}")
    print(f"    File1 value: {stats['value1_at_max_abs']:.6e}")
    print(f"    File2 value: {stats['value2_at_max_abs']:.6e}")
    print(f"    Absolute difference: {stats['max_absolute_diff']:.6e}")
    print()
    print(f"  Cell with max relative difference: {stats['cell_max_rel_diff']}")
    print(f"    File1 value: {stats['value1_at_max_rel']:.6e}")
    print(f"    File2 value: {stats['value2_at_max_rel']:.6e}")
    print(f"    Relative difference: {stats['max_relative_diff']:.6e}")
    print("=" * 60)

def main():
    if len(sys.argv) != 3:
        print("Usage: python cfd_compare.py <file1> <file2>")
        print("Example: python cfd_compare.py cpu_results.dat gpu_results.dat")
        return
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    # Compare files
    results = compare_cfd_files(file1, file2)
    
    if results is None:
        return
    
    # Check if files are identical
    if results[0] == 'identical':
        return
    
    stats, differences, relative_diff, cell_ids, values1, values2 = results
    
    # Print results
    print_comparison_results(stats)
    
    # Optionally save top 10 largest differences
    print("\nTop 10 largest absolute differences:")
    sorted_indices = np.argsort(differences)[-10:][::-1]
    for i, idx in enumerate(sorted_indices):
        print(f"{i+1:2d}. Cell {cell_ids[idx]:6d}: "
              f"Diff = {differences[idx]:.6e} "
              f"({values1[idx]:.6e} vs {values2[idx]:.6e})")

if __name__ == "__main__":
    main()
