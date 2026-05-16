# Hi-TOM Local Genotyping Pipeline

A Python reimplementation of the **Hi-TOM** pipeline (Liu et al., *Sci China Life Sci* 2019) for high-throughput CRISPR/Cas9 mutation profiling. Processes pooled amplicon reads from a 96-well plate: demultiplexes by inline barcode, extracts target windows from BWA-aligned BAMs, classifies mutations with `difflib`, and outputs per-well genotype tables in the original Hi-TOM TSV format.

Two pipeline variants are provided:

| Version | File | Description |
|---------|------|-------------|
| **v1** (recommended) | `hitom_analyze.py` | `MAPQ >= 20`, per-base Q25 on window, singleton recovery, artifact scrub |
| **v2** (original paper thresholds) | `hitom_analyze_v2.py` | `MAPQ >= 30`, barcode Q30 gate, `>10% N` / `>50% Q<=5` pre-filter, no singleton recovery, no artifact scrub |

---

## Pipeline Steps

1. **Demultiplex** — reads FASTQ headers, extracts 4-nt inline barcodes (R1 positions 5–8 for column, R2 positions 5–8 for row), fuzzy Hamming-distance-1 matching
2. **Auto-detect target window** — scans BAM coverage profile, locates peak, expands left/right until coverage drops to ~5 % of peak
3. **Extract target window** — walks CIGAR strings, extracts bases overlapping `[start, end]`
4. **Call variants** — `difflib.SequenceMatcher` with free-end-gap treatment + Hamming-distance SNP fallback
5. **Collapse haplotypes** — 5 % per-sample abundance threshold, genotype classification (`AA`, `Aa`, `aa`, `chimeric`, `-`)
6. **Output** — per-well sequence table (`*-sequence.tsv`) and 8x12 genotype matrix (`*-genotype.tsv`)

---

## Requirements

| Dependency | Version | Notes |
|-----------|---------|-------|
| Python | >= 3.8 | Standard library only (`gzip`, `difflib`, `re`, `subprocess`, `itertools`, `collections`) |
| WSL | any | Required on Windows for `samtools` and `bwa` |
| `samtools` | >= 1.10 | `samtools view`, `samtools depth`, `samtools sort`, `samtools index` |
| `bwa` | >= 0.7.17 | Must be run beforehand to generate sorted & indexed BAM |

---

## Input Files

For each dataset:

1. **`*_raw_1.fq.gz`** — R1 FASTQ
2. **`*_raw_2.fq.gz`** — R2 FASTQ
3. **`*_sorted.bam`** + **`*.bai`** — BWA-mem aligned, sorted, indexed BAM
4. **`*_reference_sequence.fasta`** — single-entry reference FASTA

---

## Usage

### 1. Align reads

```bash
bwa index <ref>.fasta
bwa mem <ref>.fasta <r1>.fq.gz <r2>.fq.gz | samtools sort -o <dataset>_sorted.bam
samtools index <dataset>_sorted.bam
```

### 2. Run the pipeline

```python
from hitom_analyze import process_dataset

process_dataset(
    dataset_num=1,
    r1_path='01_raw_1.fq.gz',
    r2_path='01_raw_2.fq.gz',
    bam_path='01_sorted.bam',
    ref_path='01_reference_sequence.fasta',
    out_seq='01-sequence.tsv',
    out_geno='01-genotype.tsv'
)
```

### 3. Validate against reference examples

```bash
python compare_final.py
```

---

## Quick Start with Demo Data

The `demo/` folder contains reduced datasets (10 % random subsample of read pairs) for quick testing:

```
demo/
  01_raw_1.fq.gz   01_reference_sequence.fasta   01_sorted.bam   01_sorted.bam.bai
  01_raw_2.fq.gz   01-sequence.tsv               01-genotype.tsv
  27_raw_1.fq.gz   27_reference_sequence.fasta   27_sorted.bam   27_sorted.bam.bai
  27_raw_2.fq.gz   27-sequence.tsv               27-genotype.tsv
```

These have already been aligned and processed. To re-run:

```python
from hitom_analyze import process_dataset
process_dataset(1, 'demo/01_raw_1.fq.gz', 'demo/01_raw_2.fq.gz',
                'demo/01_sorted.bam', 'demo/01_reference_sequence.fasta',
                'demo/01-sequence.tsv', 'demo/01-genotype.tsv')
```

---

## Output Format

### Sequence table (`*-sequence.tsv`)

| Column | Description |
|--------|-------------|
| Sort | Rank by read count |
| Reads number | Haplotype abundance |
| Ratio | Percentage of total |
| Left/Right variation type | `WT`, `nI`, `nD`, `SNP`, `Large_Indel` |
| Left/Right variation | Detail string (e.g. `A`, `GA`, `T->G`, `----`) |
| Left/Right reads seq | Extracted target-window sequence |

### Genotype table (`*-genotype.tsv`)

8 rows (A–H) x 12 columns (1–12) matrix. Genotypes:

| Code | Meaning |
|------|---------|
| `AA` | Wild-type |
| `Aa` | Heterozygous |
| `aa` | Homozygous mutant |
| `chimeric` | 3+ variant types |
| `-` | Missing data |
| `*` | In-frame indel |
| `#` | SNP |

---

## v1 vs v2 Filtering Differences

The original Hi-TOM web platform (Perl 5.16.3, BWA-MEM 0.7.10) was never open-sourced. The paper describes these filtering criteria:

| Criterion | Original Perl (paper) | v1 (`hitom_analyze.py`) | v2 (`hitom_analyze_v2.py`) |
|---|---|---|---|
| Barcode quality | Each base Phred > 30 | Hamming-distance only | Each base Phred > 30 |
| Read pre-filter | >10% Ns or >50% bases Q<=5 | None | >10% Ns or >50% bases Q<=5 |
| MAPQ threshold | >= 30 | >= 20 | >= 30 |
| Per-base Q25 on window | None | Enforced (rejects read) | None |
| Singleton recovery | No | Yes | No |
| Artifact scrub | No | Yes | No |

Neither v1 nor v2 exactly replicates the original because:
- The original variant caller is a proprietary Perl script (not `difflib`)
- The original uses BWA-MEM **0.7.10** specifically (different alignment behaviour)
- The original's duplicate-removal algorithm for amplicon data is unspecified

v1 achieves higher accuracy on the example datasets and is therefore recommended.

---

## Validation Results

Pipeline validated against the original paper's two example datasets (full, non-demo data):

| Dataset | Gene | Genotype Match | Notes |
|---------|------|---------------|-------|
| **01** | Gene1 | **100.0 %** (96/96) | Perfect agreement |
| **27** | Gene2 | **90.6 %** (87/96) | Borderline threshold cases near 5 % abundance |

Dataset 27 mismatches are exclusively in borderline cases where variant abundance hovers near the 5 % threshold. The differences stem from the original's different variant caller and quality-filtering strategy — they are not pipeline bugs.

---

## File Structure

```
.
  .gitignore                   # Excludes test-files/ (large raw data)
  hitom_analyze.py             # v1 pipeline (recommended)
  hitom_analyze_v2.py          # v2 pipeline (original paper thresholds)
  hitom_gui.py                 # Optional GUI wrapper
  compare_final.py             # Validation tool
  README.md                    # This file
  Design_info.md               # Barcode & primer definitions
  Hi-Tom/                      # Original paper & figures
  example_results/             # Reference example TSV outputs
  demo/                        # 10 % subsampled demo dataset
  test-files/                  # Full raw data (gitignored)
```

---

## Citation

> **Liu Q, et al.** "Hi-TOM: a platform for high-throughput tracking of mutations induced by CRISPR/Cas systems." *Science China Life Sciences*, 2019. doi:10.1007/s11427-018-9402-9

This is an independent reimplementation. It is not affiliated with the original authors or the Hi-TOM web platform.
