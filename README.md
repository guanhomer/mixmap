# mixmap  
**MIXture Methylation Analysis in Plasma**

mixmap provides tools for extracting and analyzing fragment-level CpG methylation patterns from Enzymatic Methyl-seq (EM-seq) or bisulfite sequencing data. It is designed for applications involving heterogeneous DNA mixtures, such as circulating cell-free DNA (cfDNA), where per-fragment methylation information is required.

The core functionality converts aligned reads (BAM) into per-fragment CpG methylation calls, enabling downstream analyses of methylation heterogeneity and mixture composition.

---

## Features

- Per-read, per-CpG methylation calling
- Fragment-level collapse of overlapping CpG observations
- Strand-aware interpretation compatible with common bisulfite aligners
- External sorting for scalability to large BAM files
- Minimal dependencies and streaming-friendly design

---

## Requirements

- Python ≥ 3.8  
- pysam  
- Unix environment with `sort` (GNU coreutils recommended)

---

## Installation

```bash
git clone https://github.com/<your-username>/mixmap.git
cd mixmap
pip install pysam
```

---

## Usage

```bash
python per_read_per_cpg.py reference.fa input.bam output.tsv
```

### Required arguments

- `reference.fa`  
  Reference genome in FASTA format (must be indexed with `.fai`)

- `input.bam`  
  Input BAM file (use `-` for stdin)

- `output.tsv`  
  Output file with per-fragment CpG methylation calls

---

## Options

```bash
-p, --min-phred INT        Minimum base quality (default: 5)
-q, --min-mapq INT         Minimum mapping quality (default: 10)

--require-flags INT        Require BAM flags (decimal or hex)
--ignore-flags INT         Ignore BAM flags (default: 0x904)

--overlap-policy {agree,best}
                           Resolve conflicting calls within fragment:
                             agree → keep only if consistent (default)
                             best  → keep highest-quality observation

--tempdir PATH             Temporary directory for intermediate files
```

---

## Output format

| Column           | Description |
|------------------|-------------|
| fragment_id      | Sequential fragment index |
| chrom            | Chromosome |
| ref_pos          | 1-based CpG coordinate |
| strand           | Genomic strand (`+` or `-`) |
| call             | Methylation state (`m` or `u`) |
| fragment_span    | Fragment length in reference coordinates |

---

## Method overview

1. Extraction (Stage 1): identify CpG-overlapping positions and generate per-read methylation calls  
2. Sorting: external sort by fragment and genomic position  
3. Collapse (Stage 2): resolve overlapping calls within fragments  

---

## Notes

- Chromosome filtering is limited to canonical human chromosomes (1–22, X, Y, MT).
- Designed for bisulfite sequencing alignments (e.g. Bismark, bwa-meth).
- Memory usage is controlled via external sorting.
