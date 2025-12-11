# HDF5 Compression Tests for EmpCylSL

## Overview

This directory contains tests for the HDF5 compression functionality implemented in the EmpCylSL cylindrical disk basis cache system.

## Test File: `cyl_compression_test.py`

This comprehensive test suite validates the HDF5 compression implementation by testing:

### Test 1: Compression Disabled
- Verifies that compression can be disabled by setting `compress: 0` in the YAML configuration
- Creates a basis with no compression and verifies the cache file is created and readable

### Test 2: Compression Enabled
- Tests multiple compression levels (1, 5, 9)
- Verifies that each compression level works correctly
- Ensures cache files are created and readable for all compression levels

### Test 3: Compressed File Readback
- Creates a compressed cache file with compression level 5
- Verifies the file can be read back correctly by a second basis instance
- Ensures data integrity is maintained through compression/decompression

### Test 4: Compression Level Validation
- Tests valid compression levels (0, 5, 9)
- Tests out-of-range compression levels (10, 15, 100)
- Verifies that out-of-range values are clamped to 9 as implemented in `setH5Params`

### Test 5: Shuffle Parameter
- Verifies the shuffle filter works in conjunction with compression
- Note: The shuffle parameter is enabled by default in `setH5Params` when compression is enabled
- Tests indirectly by verifying successful cache creation with compression

### Test 6: Compression Comparison
- Compares file sizes between uncompressed (level 0) and compressed (level 9) cache files
- Provides visual feedback on compression effectiveness
- Note: Small test datasets may not compress significantly

## Running the Tests

### Via CMake/CTest

After building the project:

```bash
cd build
ctest -R pyexpCylCompressionTest --output-on-failure
```

Or as part of the full test suite:

```bash
ctest -L long --output-on-failure
```

### Direct Execution

You can also run the test directly if pyEXP is built and in your PYTHONPATH:

```bash
cd tests/Disk
python3 cyl_compression_test.py
```

## Expected Output

When all tests pass, you should see output similar to:

```
======================================================================
EmpCylSL HDF5 Compression Test Suite
======================================================================

Test 1: Testing compression disabled (compress=0)...
  ✓ Compression disabled test passed

Test 2: Testing compression enabled with different levels...
  Testing compression level 1...
    ✓ Compression level 1 test passed
  Testing compression level 5...
    ✓ Compression level 5 test passed
  Testing compression level 9...
    ✓ Compression level 9 test passed
  ✓ All compression level tests passed

[... more tests ...]

======================================================================
Test Summary
======================================================================
Passed: 6/6
✓ All tests passed!
```

## Implementation Details

The compression functionality is implemented in:
- `src/Cylinder.cc`: Reads the `compress` parameter from YAML config and calls `setH5Params`
- `include/EmpCylSL.H`: Defines `setH5Params` method with validation
- `exputil/EmpCylSL.cc`: Applies compression settings when writing HDF5 cache files

The `setH5Params` method signature:
```cpp
void setH5Params(unsigned compress, bool shuffle=true, bool szip=false)
```

Where:
- `compress`: Gzip compression level 0-9 (0=off, 9=max)
- `shuffle`: Enable byte shuffle filter (default: true)
- `szip`: Enable szip compression (default: false)

## Test Coverage

These tests address the feedback from PR #183:
- ✅ Compression can be enabled/disabled via the compress parameter
- ✅ Compressed cache files can be read back correctly
- ✅ The shuffle parameter works as expected
- ✅ Validation of the compression level (0-9 range)
