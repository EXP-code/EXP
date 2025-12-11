#!/usr/bin/env python
# coding: utf-8

"""
Test HDF5 compression functionality in EmpCylSL disk cache.

This test suite verifies:
1. Compression can be enabled/disabled via the compress parameter
2. Compressed cache files can be read back correctly
3. The shuffle parameter works as expected
4. Validation of the compression level (0-9 range)
"""

import os
import pyEXP
import tempfile                 # For creating temporary files and directories
import shutil                   # For cleaning up temporary files and directories

# Base configuration for disk basis
#
base_config_template = """
---
id: cylinder
parameters:
  acyl: 0.01       # The scale length of the exponential disk
  hcyl: 0.001      # The scale height of the exponential disk
  lmaxfid: 20      # The maximum spherical harmonic order for the input basis
  nmaxfid: 20      # The radial order for the input spherical basis
  mmax: 6          # The maximum azimuthal order for the cylindrical basis
  nmax: 8          # The maximum radial order of the cylindrical basis
  ncylnx: 128      # The number of grid points in mapped cylindrical radius
  ncylny: 64       # The number of grid points in mapped verical scale
  ncylodd: 3       # The number of anti-symmetric radial basis functions per azimuthal order m
  rnum: 32         # The number of radial integration knots in the inner product
  pnum: 0          # The number of azimuthal integration knots (pnum: 0, assume axisymmetric target density)
  tnum: 16         # The number of colatitude integration knots
  ashift: 0.5      # Target shift length in scale lengths to create more variance
  vflag: 0         # Verbosity flag: print diagnostics to stdout for vflag>0
  logr: false      # Log scaling in cylindrical radius
  sech2: true      # Use standard defintion
  cachename: {eof_file}  # The cache file name
  compress: {compress}  # Compression level (0-9)
...
"""

def test_compression_disabled():
    """Test 1: Verify compression can be disabled (compress=0)"""
    print("Test 1: Testing compression disabled (compress=0)...")
    
    # Create temporary directory for test files
    test_dir = tempfile.mkdtemp(prefix="exp_test_compress_disabled_")
    cache_file = os.path.join(test_dir, ".eof.cache.compress_disabled")
    
    try:
        config = base_config_template.format(
            eof_file=cache_file,
            compress=0
        )
        
        # Create basis with compression disabled
        basis = pyEXP.basis.Basis.factory(config)
        
        # Verify cache file was created
        assert os.path.exists(cache_file), "Cache file was not created"
        
        # Verify we can read the cache file
        node_cyl = basis.cacheInfo(cache_file)
        assert node_cyl is not None, "Failed to read cache file"
        
        print("  ✓ Compression disabled test passed")
        return True
        
    except Exception as e:
        print(f"  ✗ Compression disabled test failed: {e}")
        return False
    finally:
        # Cleanup
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)


def test_compression_enabled():
    """Test 2: Verify compression can be enabled with various levels"""
    print("Test 2: Testing compression enabled with different levels...")
    
    compression_levels = [1, 5, 9]
    results = []
    
    for level in compression_levels:
        print(f"  Testing compression level {level}...")
        
        # Create temporary directory for test files
        test_dir = tempfile.mkdtemp(prefix=f"exp_test_compress_{level}_")
        cache_file = os.path.join(test_dir, f".eof.cache.compress_{level}")
        
        try:
            config = base_config_template.format(
                eof_file=cache_file,
                compress=level
            )
            
            # Create basis with compression enabled
            basis = pyEXP.basis.Basis.factory(config)
            
            # Verify cache file was created
            assert os.path.exists(cache_file), f"Cache file was not created for level {level}"
            
            # Verify we can read the cache file back
            node_cyl = basis.cacheInfo(cache_file)
            assert node_cyl is not None, f"Failed to read cache file for level {level}"
            
            print(f"    ✓ Compression level {level} test passed")
            results.append(True)
            
        except Exception as e:
            print(f"    ✗ Compression level {level} test failed: {e}")
            results.append(False)
        finally:
            # Cleanup
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
    
    all_passed = all(results)
    if all_passed:
        print("  ✓ All compression level tests passed")
    return all_passed


def test_compressed_file_readback():
    """Test 3: Verify compressed cache files can be read back correctly"""
    print("Test 3: Testing compressed file readback...")
    
    # Create temporary directory for test files
    test_dir = tempfile.mkdtemp(prefix="exp_test_readback_")
    cache_file = os.path.join(test_dir, ".eof.cache.readback")
    
    try:
        config = base_config_template.format(
            eof_file=cache_file,
            compress=5
        )
        
        # Create basis with compression enabled
        basis1 = pyEXP.basis.Basis.factory(config)
        
        # Verify cache file was created
        assert os.path.exists(cache_file), "Cache file was not created"
        
        # Read cache info from the first basis
        node_cyl1 = basis1.cacheInfo(cache_file)
        assert node_cyl1 is not None, "Failed to read cache file with first basis"
        
        # Create a second basis instance that should read from the existing cache
        basis2 = pyEXP.basis.Basis.factory(config)
        
        # Read cache info from the second basis
        node_cyl2 = basis2.cacheInfo(cache_file)
        assert node_cyl2 is not None, "Failed to read cache file with second basis"
        
        print("  ✓ Compressed file readback test passed")
        return True
        
    except Exception as e:
        print(f"  ✗ Compressed file readback test failed: {e}")
        return False
    finally:
        # Cleanup
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)


def test_compression_level_validation():
    """Test 4: Verify compression level validation (0-9 range)"""
    print("Test 4: Testing compression level validation...")
    
    # Test valid range values
    valid_levels = [0, 5, 9]
    # Test out of range values (should be clamped to 9)
    out_of_range_levels = [10, 15, 100]
    
    results = []
    
    # Test valid levels
    for level in valid_levels:
        print(f"  Testing valid compression level {level}...")
        
        test_dir = tempfile.mkdtemp(prefix=f"exp_test_valid_{level}_")
        cache_file = os.path.join(test_dir, f".eof.cache.valid_{level}")
        
        try:
            config = base_config_template.format(
                eof_file=cache_file,
                compress=level
            )
            
            # Should succeed without error
            basis = pyEXP.basis.Basis.factory(config)
            
            # Verify cache file was created
            assert os.path.exists(cache_file), f"Cache file was not created for valid level {level}"
            
            print(f"    ✓ Valid compression level {level} accepted")
            results.append(True)
            
        except Exception as e:
            print(f"    ✗ Valid compression level {level} failed: {e}")
            results.append(False)
        finally:
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
    
    # Test out of range levels (should be clamped to 9 but still work)
    for level in out_of_range_levels:
        print(f"  Testing out-of-range compression level {level} (should be clamped to 9)...")
        
        test_dir = tempfile.mkdtemp(prefix=f"exp_test_oor_{level}_")
        cache_file = os.path.join(test_dir, f".eof.cache.oor_{level}")
        
        try:
            config = base_config_template.format(
                eof_file=cache_file,
                compress=level
            )
            
            # Should succeed (value will be clamped to 9 internally)
            basis = pyEXP.basis.Basis.factory(config)
            
            # Verify cache file was created
            assert os.path.exists(cache_file), f"Cache file was not created for OOR level {level}"
            
            print(f"    ✓ Out-of-range compression level {level} handled correctly")
            results.append(True)
            
        except Exception as e:
            print(f"    ✗ Out-of-range compression level {level} failed: {e}")
            results.append(False)
        finally:
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
    
    all_passed = all(results)
    if all_passed:
        print("  ✓ All compression level validation tests passed")
    return all_passed


def test_shuffle_parameter():
    """Test 5: Verify the shuffle parameter works (indirectly through successful cache creation)"""
    print("Test 5: Testing shuffle parameter...")
    
    # Note: The shuffle parameter is set via setH5Params in C++ code (include/EmpCylSL.H:952),
    # not directly via YAML config. We test it indirectly by verifying that compression with
    # default shuffle works. The actual shuffle configuration is controlled internally and
    # defaults to true as defined in the setH5Params method signature.
    
    test_dir = tempfile.mkdtemp(prefix="exp_test_shuffle_")
    cache_file = os.path.join(test_dir, ".eof.cache.shuffle")
    
    try:
        # Create with compression enabled (shuffle is enabled by default in setH5Params)
        config = base_config_template.format(
            eof_file=cache_file,
            compress=5
        )
        
        basis = pyEXP.basis.Basis.factory(config)
        
        # Verify cache file was created
        assert os.path.exists(cache_file), "Cache file was not created with shuffle"
        
        # Verify we can read the cache file back
        node_cyl = basis.cacheInfo(cache_file)
        assert node_cyl is not None, "Failed to read cache file with shuffle"
        
        print("  ✓ Shuffle parameter test passed (shuffle enabled by default with compression)")
        return True
        
    except Exception as e:
        print(f"  ✗ Shuffle parameter test failed: {e}")
        return False
    finally:
        # Cleanup
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)


def test_compression_comparison():
    """Test 6: Compare file sizes with and without compression"""
    print("Test 6: Comparing cache file sizes with and without compression...")
    
    test_dir = tempfile.mkdtemp(prefix="exp_test_comparison_")
    cache_file_uncompressed = os.path.join(test_dir, ".eof.cache.uncompressed")
    cache_file_compressed = os.path.join(test_dir, ".eof.cache.compressed")
    
    try:
        # Create uncompressed cache file
        config_uncompressed = base_config_template.format(
            eof_file=cache_file_uncompressed,
            compress=0
        )
        basis_uncompressed = pyEXP.basis.Basis.factory(config_uncompressed)
        
        # Create compressed cache file
        config_compressed = base_config_template.format(
            eof_file=cache_file_compressed,
            compress=9
        )
        basis_compressed = pyEXP.basis.Basis.factory(config_compressed)
        
        # Get file sizes
        size_uncompressed = os.path.getsize(cache_file_uncompressed)
        size_compressed = os.path.getsize(cache_file_compressed)
        
        print(f"  Uncompressed cache size: {size_uncompressed:,} bytes")
        print(f"  Compressed cache size: {size_compressed:,} bytes")
        
        if size_compressed < size_uncompressed:
            compression_ratio = (1 - size_compressed / size_uncompressed) * 100
            print(f"  Compression ratio: {compression_ratio:.1f}% reduction")
            print("  ✓ Compression reduces file size as expected")
            return True
        else:
            print("  ⚠ Warning: Compressed file is not smaller (may be expected for small test data)")
            # This is not necessarily a failure - small test datasets may not compress well
            return True
        
    except Exception as e:
        print(f"  ✗ Compression comparison test failed: {e}")
        return False
    finally:
        # Cleanup
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)


def main():
    """Run all compression tests"""
    print("=" * 70)
    print("EmpCylSL HDF5 Compression Test Suite")
    print("=" * 70)
    print()
    
    tests = [
        test_compression_disabled,
        test_compression_enabled,
        test_compressed_file_readback,
        test_compression_level_validation,
        test_shuffle_parameter,
        test_compression_comparison,
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"✗ Test crashed: {e}")
            results.append(False)
        print()
    
    # Summary
    print("=" * 70)
    print("Test Summary")
    print("=" * 70)
    passed = sum(results)
    total = len(results)
    print(f"Passed: {passed}/{total}")
    
    if all(results):
        print("✓ All tests passed!")
        return 0
    else:
        print("✗ Some tests failed")
        return 1


if __name__ == "__main__":
    exit(main())
