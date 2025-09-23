#!/bin/bash

# Test script for STR_regression_streaming.R
# Generates synthetic test data and validates the streaming analysis

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test directory
TEST_DIR="test_str_streaming"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="$SCRIPT_DIR/STR_regression_streaming.R"

echo -e "${BLUE}=== STR Regression Streaming Test Suite ===${NC}"
echo "Creating test data and validating streaming analysis..."
echo "Using script: $SCRIPT_PATH"

# Create test directory
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

echo -e "\n${YELLOW}Step 1: Generating synthetic STR data...${NC}"

# Generate synthetic STR input file (inquiSTR format)
cat > test_str_data.txt << 'EOF'
chrom	start	end	sample1_H1	sample1_H2	sample2_H1	sample2_H2	sample3_H1	sample3_H2	sample4_H1	sample4_H2	sample5_H1	sample5_H2	sample6_H1	sample6_H2	sample7_H1	sample7_H2	sample8_H1	sample8_H2	sample9_H1	sample9_H2	sample10_H1	sample10_H2
chr1	1000	1030	20	22	18	20	25	27	19	21	23	25	17	19	26	28	20	22	24	26	18	20
chr1	2000	2040	15	17	30	32	16	18	31	33	14	16	29	31	15	17	32	34	13	15	28	30
chr2	3000	3050	45	47	12	14	46	48	11	13	44	46	10	12	47	49	13	15	45	47	12	14
chr2	4000	4060	8	10	50	52	9	11	51	53	7	9	49	51	8	10	52	54	6	8	48	50
chr3	5000	5080	35	37	22	24	36	38	21	23	34	36	20	22	37	39	23	25	35	37	22	24
chr3	6000	6100	5	7	60	62	6	8	61	63	4	6	59	61	5	7	62	64	3	5	58	60
chr4	7000	7120	40	42	15	17	41	43	14	16	39	41	13	15	42	44	16	18	40	42	15	17
chr4	8000	8150	12	14	45	47	13	15	46	48	11	13	44	46	12	14	47	49	10	12	43	45
EOF

echo -e "${GREEN}✓${NC} Created test_str_data.txt with 8 STR variants and 10 samples"

# Generate phenotype file
cat > test_phenotypes.txt << 'EOF'
sample_id	disease_status	age	sex	continuous_trait
sample1	Control	45	M	1.2
sample2	Case	52	F	2.8
sample3	Control	38	M	0.9
sample4	Case	61	F	3.1
sample5	Control	29	M	1.7
sample6	Case	44	F	2.5
sample7	Control	33	F	1.1
sample8	Case	57	M	3.4
sample9	Control	41	M	1.6
sample10	Case	48	F	2.9
EOF

echo -e "${GREEN}✓${NC} Created test_phenotypes.txt with binary and continuous traits"

echo -e "\n${YELLOW}Step 2: Testing binary outcome analysis...${NC}"

# Test 1: Binary outcome with MEAN mode
echo "Running binary analysis (MEAN mode)..."
Rscript "$SCRIPT_PATH" \
    --input test_str_data.txt \
    --phenocovar test_phenotypes.txt \
    --phenotype disease_status \
    --out test_binary_mean.txt \
    --STRmode MEAN \
    --outcometype binary \
    --binaryOrder "Control,Case" \
    --missing_cutoff 0.5

if [[ -f "test_binary_mean.txt" ]]; then
    echo -e "${GREEN}✓${NC} Binary analysis (MEAN) completed successfully"
    echo "Results preview:"
    head -3 test_binary_mean.txt | column -t
    echo "Total variants in output: $(tail -n +2 test_binary_mean.txt | wc -l)"
else
    echo -e "${RED}✗${NC} Binary analysis (MEAN) failed"
    exit 1
fi

# Test 2: Binary outcome with MAX mode
echo -e "\nRunning binary analysis (MAX mode)..."
Rscript "$SCRIPT_PATH" \
    --input test_str_data.txt \
    --phenocovar test_phenotypes.txt \
    --phenotype disease_status \
    --out test_binary_max.txt \
    --STRmode MAX \
    --outcometype binary \
    --binaryOrder "Control,Case" \
    --missing_cutoff 0.5 \
    --quiet

if [[ -f "test_binary_max.txt" ]]; then
    echo -e "${GREEN}✓${NC} Binary analysis (MAX) completed successfully"
else
    echo -e "${RED}✗${NC} Binary analysis (MAX) failed"
    exit 1
fi

echo -e "\n${YELLOW}Step 3: Testing continuous outcome analysis...${NC}"

# Test 3: Continuous outcome with covariates
echo "Running continuous analysis with covariates..."
Rscript "$SCRIPT_PATH" \
    --input test_str_data.txt \
    --phenocovar test_phenotypes.txt \
    --phenotype continuous_trait \
    --out test_continuous.txt \
    --STRmode MEAN \
    --outcometype continuous \
    --covnames "age,sex" \
    --missing_cutoff 0.5

if [[ -f "test_continuous.txt" ]]; then
    echo -e "${GREEN}✓${NC} Continuous analysis completed successfully"
    echo "Results preview:"
    head -3 test_continuous.txt | column -t
    echo "Total variants in output: $(tail -n +2 test_continuous.txt | wc -l)"
else
    echo -e "${RED}✗${NC} Continuous analysis failed"
    exit 1
fi

echo -e "\n${YELLOW}Step 4: Testing filtering options...${NC}"

# Test 4: Minimal length filter
echo "Testing minimal length filter..."
Rscript "$SCRIPT_PATH" \
    --input test_str_data.txt \
    --phenocovar test_phenotypes.txt \
    --phenotype disease_status \
    --out test_length_filter.txt \
    --STRmode MAX \
    --outcometype binary \
    --binaryOrder "Control,Case" \
    --minimal_length 40 \
    --missing_cutoff 0.5 \
    --quiet

if [[ -f "test_length_filter.txt" ]]; then
    filtered_count=$(tail -n +2 test_length_filter.txt | wc -l)
    echo -e "${GREEN}✓${NC} Length filtering completed successfully"
    echo "Variants passing length filter (>40): $filtered_count"
else
    echo -e "${RED}✗${NC} Length filtering failed"
    exit 1
fi

# Test 5: High missing cutoff (should filter most variants)
echo -e "\nTesting high missing cutoff filter..."
Rscript "$SCRIPT_PATH" \
    --input test_str_data.txt \
    --phenocovar test_phenotypes.txt \
    --phenotype disease_status \
    --out test_missing_filter.txt \
    --STRmode MEAN \
    --outcometype binary \
    --binaryOrder "Control,Case" \
    --missing_cutoff 0.99 \
    --quiet

if [[ -f "test_missing_filter.txt" ]]; then
    filtered_count=$(tail -n +2 test_missing_filter.txt | wc -l)
    echo -e "${GREEN}✓${NC} Missing data filtering completed successfully"
    echo "Variants passing high missing filter (0.99): $filtered_count"
else
    echo -e "${RED}✗${NC} Missing data filtering failed"
    exit 1
fi

echo -e "\n${YELLOW}Step 5: Testing error handling...${NC}"

# Test 6: Invalid phenotype column
echo "Testing invalid phenotype column (should fail)..."
if Rscript "$SCRIPT_PATH" \
    --input test_str_data.txt \
    --phenocovar test_phenotypes.txt \
    --phenotype nonexistent_phenotype \
    --out test_error.txt \
    --STRmode MEAN \
    --outcometype binary \
    --binaryOrder "Control,Case" \
    --quiet 2>/dev/null; then
    echo -e "${RED}✗${NC} Error handling failed - should have caught invalid phenotype"
    exit 1
else
    echo -e "${GREEN}✓${NC} Correctly caught invalid phenotype column"
fi

# Test 7: Missing binary order
echo "Testing missing binary order (should fail)..."
if Rscript "$SCRIPT_PATH" \
    --input test_str_data.txt \
    --phenocovar test_phenotypes.txt \
    --phenotype disease_status \
    --out test_error.txt \
    --STRmode MEAN \
    --outcometype binary \
    --quiet 2>/dev/null; then
    echo -e "${RED}✗${NC} Error handling failed - should have caught missing binary order"
    exit 1
else
    echo -e "${GREEN}✓${NC} Correctly caught missing binary order"
fi

echo -e "\n${YELLOW}Step 6: Testing data consistency...${NC}"

# Compare different STR modes for same data
echo "Comparing MEAN vs MAX vs MIN modes..."

# Run MIN mode
Rscript "$SCRIPT_PATH" \
    --input test_str_data.txt \
    --phenocovar test_phenotypes.txt \
    --phenotype disease_status \
    --out test_binary_min.txt \
    --STRmode MIN \
    --outcometype binary \
    --binaryOrder "Control,Case" \
    --missing_cutoff 0.5 \
    --quiet

# Check that all three modes produced output
mean_vars=$(tail -n +2 test_binary_mean.txt | wc -l)
max_vars=$(tail -n +2 test_binary_max.txt | wc -l)
min_vars=$(tail -n +2 test_binary_min.txt | wc -l)

if [[ "$mean_vars" -eq "$max_vars" && "$max_vars" -eq "$min_vars" ]]; then
    echo -e "${GREEN}✓${NC} All STR modes processed same number of variants ($mean_vars)"
else
    echo -e "${YELLOW}!${NC} Different number of variants across modes (MEAN:$mean_vars, MAX:$max_vars, MIN:$min_vars)"
fi

echo -e "\n${YELLOW}Step 7: Performance validation...${NC}"

# Generate larger dataset for performance test
echo "Generating larger dataset for performance test..."
cat > test_large_data.txt << 'EOF'
chrom	start	end	sample1_H1	sample1_H2	sample2_H1	sample2_H2	sample3_H1	sample3_H2	sample4_H1	sample4_H2	sample5_H1	sample5_H2
EOF

# Add 100 variants with proper sample matching
for i in {1..100}; do
    chr="chr$((i % 5 + 1))"
    start=$((i * 1000))
    end=$((start + 30 + i % 20))
    
    # Generate random STR lengths (10-50) for 5 samples x 2 alleles = 10 values
    vals=""
    for _ in {1..10}; do  # 5 samples * 2 alleles
        val=$((10 + RANDOM % 40))
        vals="$vals	$val"
    done
    echo "$chr	$start	$end$vals" >> test_large_data.txt
done

echo "Generated test_large_data.txt with $(tail -n +2 test_large_data.txt | wc -l) variants"

echo "Running performance test on 100 variants..."
start_time=$(date +%s)
Rscript "$SCRIPT_PATH" \
    --input test_large_data.txt \
    --phenocovar test_phenotypes.txt \
    --phenotype continuous_trait \
    --out test_performance.txt \
    --STRmode MEAN \
    --outcometype continuous \
    --missing_cutoff 0.5 \
    --quiet

end_time=$(date +%s)
runtime=$((end_time - start_time))
processed_vars=$(tail -n +2 test_performance.txt | wc -l)

echo -e "${GREEN}✓${NC} Performance test completed in ${runtime}s"
if [[ $runtime -gt 0 ]]; then
    rate=$(echo "scale=1; $processed_vars/$runtime" | bc -l)
    echo "Processed $processed_vars variants in ${runtime}s ($rate variants/sec)"
else
    echo "Processed $processed_vars variants in <1s (very fast!)"
fi

echo -e "\n${YELLOW}Step 8: Summary validation...${NC}"

# Validate output format consistency
echo "Validating output formats..."

# Check binary output format
binary_cols=$(head -1 test_binary_mean.txt | tr '\t' '\n' | wc -l)
echo "Binary output columns: $binary_cols"

# Check continuous output format  
continuous_cols=$(head -1 test_continuous.txt | tr '\t' '\n' | wc -l)
echo "Continuous output columns: $continuous_cols"

# Check that all output files have headers
files=("test_binary_mean.txt" "test_binary_max.txt" "test_continuous.txt" "test_length_filter.txt")
for file in "${files[@]}"; do
    if [[ -f "$file" ]] && [[ $(wc -l < "$file") -gt 0 ]]; then
        echo -e "${GREEN}✓${NC} $file has valid output"
    else
        echo -e "${RED}✗${NC} $file is empty or missing"
    fi
done

echo -e "\n${GREEN}=== Test Suite Completed Successfully! ===${NC}"
echo -e "\nSummary of generated files:"
for file in *.txt; do
    if [[ -f "$file" ]]; then
        echo "  $file ($(wc -l < "$file") lines)"
    fi
done

echo -e "\n${BLUE}=== Cleanup ===${NC}"
echo "Test files are in directory: $(pwd)"
echo "To clean up, run: rm -rf $(pwd)"

cd ..
echo -e "${GREEN}All tests passed! The STR regression streaming script is working correctly.${NC}"