# STR Regression Streaming Test Suite

This directory contains a comprehensive test suite for the streaming STR regression analysis script.

## Files

- `test_str_streaming.sh` - Main test script that generates synthetic data and validates the streaming analysis
- `STR_regression_streaming.R` - The streaming STR regression analysis script being tested

## What the test does

The test suite performs the following validations:

### 1. **Data Generation**

- Creates synthetic STR data with 8 variants across 10 samples
- Generates phenotype data with binary and continuous traits
- Uses realistic STR length distributions

### 2. **Core Functionality Tests**

- ✅ Binary outcome analysis (MEAN/MAX/MIN modes)
- ✅ Continuous outcome analysis with covariates
- ✅ Missing data filtering
- ✅ Minimal length filtering
- ✅ Different STR combination modes

### 3. **Error Handling Tests**

- ✅ Invalid phenotype column names
- ✅ Missing required parameters
- ✅ File existence validation

### 4. **Performance Validation**

- ✅ Processing 100 variants for performance metrics
- ✅ Memory usage tracking
- ✅ Output format consistency

### 5. **Data Consistency**

- ✅ Compares results across different STR modes
- ✅ Validates output file formats
- ✅ Checks filtering effectiveness

## Running the tests

```bash
# Run the complete test suite
./test_str_streaming.sh

# The script will create a test_str_streaming directory with all test files
```

## Expected Output

The test suite validates:

- Binary regression results with OR, confidence intervals, and p-values
- Continuous regression results with Beta coefficients
- Proper filtering based on missing data and length thresholds
- Error handling for invalid inputs
- Performance metrics for larger datasets

## Test Data Format

**STR Input Format:**

```text
chrom   start   end   sample1_H1   sample1_H2   sample2_H1   sample2_H2   ...
chr1    1000    1030  20           22           18           20           ...
```

**Phenotype Format:**

```text
sample_id    disease_status    age    sex    continuous_trait
sample1      Control          45     M      1.2
sample2      Case             52     F      2.8
```

## Cleanup

After testing, remove the test directory:

```bash
rm -rf test_str_streaming
```
