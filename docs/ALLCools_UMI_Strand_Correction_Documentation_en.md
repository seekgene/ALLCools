# ALLCools Strand Correction and UMI Correction Mechanisms Explained

## Overview

ALLCools is a toolkit for single-cell methylation data analysis that includes two important data processing mechanisms:
1. **Strand Correction** - Unifying strand orientation representation from different alignment tools (First Step)
2. **UMI Correction** - Removing PCR duplicates and correcting sequencing errors (Second Step)

This document provides detailed explanations of the principles, implementation, and applications of these two mechanisms. **Note: ALLCools processes data in a specific order: strand correction first, then UMI correction.**

## 1. Strand Correction Mechanism (First Step Processing)

### 1.1 Background Problem

Different BS-seq alignment tools handle strand orientation differently:

#### Bismark Approach
- C→T converted reads are forced to map to the positive strand (+)
- G→A converted reads are forced to map to the negative strand (-)
- Uses XG tag to mark conversion type (XG="CT" or XG="GA")

#### Hisat-3n Approach
- Maintains original orientation of R1/R2
- Same conversion type may have different strand orientations
- Uses YZ tag to mark conversion type (YZ="+" or YZ="-")

### 1.2 ALLCools Solution

ALLCools unifies strand orientation through the `_convert_bam_strandness` function:

```python
def _convert_bam_strandness(in_bam_path, out_bam_path):
    for read in in_bam:
        # Auto-detect tag type
        if read.has_tag("YZ"):
            is_ct_func = _is_read_ct_conversion_hisat3n
        elif read.has_tag("XG"):
            is_ct_func = _is_read_ct_conversion_bismark
        
        # Unify strand orientation
        if is_ct_func(read):  # C→T conversion
            read.is_forward = True
        else:  # G→A conversion
            read.is_forward = False
```

### 1.3 Tag Detection Functions

```python
def _is_read_ct_conversion_hisat3n(read):
    """Hisat-3n: YZ="+" indicates C→T conversion"""
    return read.get_tag("YZ") == "+"

def _is_read_ct_conversion_bismark(read):
    """Bismark: XG="CT" indicates C→T conversion"""
    return read.get_tag("XG") == "CT"
```

### 1.4 Impact on Downstream Analysis

After strand unification, mpileup analysis becomes consistent:
- **C sites**: Only count positive strand reads
  - `.` indicates unmethylated (C remains unchanged)
  - `T` indicates methylated (C→T conversion)
- **G sites**: Only count negative strand reads
  - `,` indicates unmethylated (G remains unchanged)  
  - `a` indicates methylated (G→A conversion)

## 2. UMI Correction Mechanism (Second Step Processing)

### 2.1 Background

In single-cell methylation sequencing, UMI (Unique Molecular Identifier) is used for:
- Identifying original molecules, distinguishing PCR duplicates
- Correcting sequencing errors through consensus checking of multiple reads
- Improving data quality and quantification accuracy

**Important**: UMI correction is performed after strand unification, ensuring consistency of mpileup data.

### 2.2 Input Data Format

UMI correction processes samtools mpileup output in the following format:
```
chr1    46930353    C    5    .T.T,    IIHII    GAAGGTGTGTAT,GAAGGTGTGTAT,GAAGGTGTGTAC,GAAGGTGTGTAT,GAAGGTGTGTAT
```

Column meanings:
- Column 1: Chromosome
- Column 2: Position
- Column 3: Reference base
- Column 4: Coverage depth
- Column 5: Base sequence (. indicates same as reference, letters indicate differences)
- Column 6: Quality scores
- Column 7: UMI sequences (comma-separated)

### 2.3 Correction Workflow

#### Step 1: Filter Relevant Bases
- **C sites**: Only keep `.` (unmethylated) and `T` (C→T conversion) reads
- **G sites**: Only keep `,` (unmethylated) and `a` (G→A conversion) reads

#### Step 2: UMI Merging (Based on Hamming Distance)
```python
def hamming_distance(str1, str2):
    """Calculate Hamming distance between two strings"""
    return sum(ch1 != ch2 for ch1, ch2 in zip(str1, str2))
```

Merging rules:
- UMIs with Hamming distance ≤ 1 are considered from the same original molecule
- When merging, keep the UMI with higher count
- When counts are equal, keep both UMIs without merging

#### Step 3: Intra-UMI Base Correction
For multiple reads within each UMI, perform base correction:

1. **Count Priority**: Select the base with the highest occurrence count
2. **Quality Score Tiebreaker**: When counts are equal, select the base with highest average quality score
3. **Unified Correction**: Unify all reads in that UMI to the corrected base

```python
# Quality score conversion
def phred33_to_quality(ascii_char):
    return ord(ascii_char) - 33

# Base correction logic
if len(tied_bases) == 1:
    corrected_base = most_common_base
else:
    # Calculate average quality scores
    for base, qualities in base_quality_map.items():
        avg_quality = sum(qualities) / len(qualities)
        if avg_quality > best_avg_quality:
            best_base = base
```

#### Step 4: Output Results
- Each UMI keeps only one representative read
- Coverage depth equals the number of unique UMIs
- All bases and quality scores are post-correction results

### 2.4 Correction Example

**Input**:
```
chr1  46930353  C  5  .T.T,  IIHII  GAAGGTGTGTAT,GAAGGTGTGTAT,GAAGGTGTGTAC,GAAGGTGTGTAT,GAAGGTGTGTAT
```

**Processing Steps**:
1. UMI counting: GAAGGTGTGTAT(4), GAAGGTGTGTAC(1)
2. Hamming distance: GAAGGTGTGTAT vs GAAGGTGTGTAC = 1
3. UMI merging: Keep GAAGGTGTGTAT (higher count)
4. Base correction: . (3 times) > T (2 times), select `.`
5. Unified correction: All bases changed to `.`

**Output**:
```
chr1  46930353  C  1  .  I  GAAGGTGTGTAT
```

## 3. Practical Applications

### 3.1 Use Cases

1. **Single-cell methylation analysis**: Requires strand unification and UMI deduplication/correction
2. **Multi-aligner compatibility**: Unifies output formats from Bismark and Hisat-3n
3. **Data quality control**: Improves data accuracy through UMI correction

### 3.2 Processing Order

Complete ALLCools processing workflow:
1. **Step 1: Strand Conversion** - Use `convert_bam_strandness=True`
2. **Step 2: UMI Correction** - Use `tag="UMI"`

### 3.3 Parameter Configuration

In the `bam_to_allc` function:
```python
bam_to_allc(
    bam_path="input.bam",
    reference_fasta="ref.fa", 
    convert_bam_strandness=True,  # Step 1: Enable strand correction
    tag="UMI"  # Step 2: Enable UMI correction
)
```

### 3.4 Output Files

- `output.allc.tsv.gz`: Main ALLC format file
- `output_mpl_old.txt`: Original mpileup output
- `output_mpl_correction.txt`: UMI-corrected mpileup output
- `output.allc.tsv.gz.count.csv`: Statistical information

## 4. Advantages and Limitations

### 4.1 Advantages

**Strand Correction**:
- Compatible with multiple alignment tools
- Unified analysis workflow
- Automatic detection and processing
- Maintains data integrity
- Provides consistent data foundation for subsequent UMI correction

**UMI Correction**:
- Effectively removes PCR duplicates
- Corrects sequencing errors
- Improves quantification accuracy
- Supports quality score weighting
- Reliable correction based on unified strand orientation

### 4.2 Limitations

**Strand Correction**:
- Depends on specific tags (XG/YZ)
- Does not support other alignment tools
- Requires additional processing steps

**UMI Correction**:
- Requires UMI information
- High computational complexity
- May over-correct low-coverage regions
- Depends on strand unification as prerequisite

## 5. Best Practices

1. **Data Preprocessing**: Ensure BAM files contain correct conversion tags
2. **Processing Order**: Strictly follow strand conversion → UMI correction sequence
3. **Parameter Tuning**: Adjust quality thresholds and coverage requirements based on data characteristics
4. **Quality Control**: Check statistical information before and after each step
5. **Result Validation**: Compare consistency of results from different processing approaches

## 6. Summary

ALLCools' strand correction and UMI correction mechanisms provide powerful data preprocessing capabilities for single-cell methylation data analysis. **The key is the correct processing order**: first perform strand unification to ensure mpileup data consistency, then perform UMI correction, which significantly improves data quality and analysis accuracy through intelligent UMI merging and base correction. This two-step processing strategy fully considers various situations in practical applications, providing researchers with reliable analysis tools.