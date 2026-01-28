# NUMT_Analysis

A bioinformatics pipeline for the detection and characterization of nuclear-embedded mitochondrial DNA segments (NUMTs) utilizing pangenome graphs.

---

## 1. Installation

This workflow uses [pixi](https://prefix.dev/) for reproducible environment management.

```bash
# Clone the repository
git clone https://github.com/LiantingFu/NUMT_Analysis.git
cd NUMT_Analysis

# Install project environments/dependencies using pixi
# (ensure pixi is installed on your system first)
pixi install
```

---

## 2. Configuration & Inputs

The pipeline requires a configuration file (`config.json`) and a metadata table for mitochondrial references (`chrM_info`).

### 2.1 Configuration (`config.json`)

Create a `config.json` file in the repository root specifying input paths. `ref` specifies the reference genome used for coordinate definition. Example:

```json
{
  "ref": "path/to/reference_genome.fasta",
  "chrM_info": "path/to/chrM_info.tsv",
  "pangenome_vcf": "path/to/pangenome_derived.vcf",
  "prefix": "output_prefix_name"
}
```

### 2.2 Mitochondrial metadata (`chrM_info`)

`chrM_info` is a TSV file describing mitochondrial references used during NUMT detection. Required columns:

| Column            | Description                                                | Example                                          |
| ----------------- | ---------------------------------------------------------- | ------------------------------------------------ |
| **species**       | Unique identifier for the sample or species                | `human_chrM_1`                                   |
| **fasta**         | Path to the raw mitochondrial FASTA file                   | `chrM_datasets/human_chrM.1.fasta`               |
| **repeat2_fasta** | Path to the mitochondrial FASTA with the sequence duplicated head-to-tail (2x) to handle circularity  | `chrM_datasets/human_chrM_1.chrM_combined.fasta` |
| **chain**         | Path to the alignment chain file (e.g., alignment to rCRS) | `chrM_datasets/human_chrM_1.aln.rCRS.chain`      |

---

## 3. Usage

Run the workflow using **Snakemake**. Adjust the number of threads (`-j`) according to your available resources. Example:

```bash
snakemake -p -j <threads> --configfile config.json
```

---

## 4. Output

Results are written to `results/final-output/`. Example output structure:

```text
results/
└── final-output
    ├── <prefix>.NUMT_detection_results.tsv         # NUMT detection report
    ├── <prefix>.NUMT_detection_results.vcf         # NUMT-associated variants (VCF)
    ├── <prefix>.fixed_reference-NUMT.info.tsv      # Details of fixed reference NUMTs
    ├── <prefix>.polymorphic_pangenome-NUMT.info.tsv # Details of polymorphic pangenome NUMTs
    └── <prefix>.polymorphic_reference-NUMT.info.tsv # Details of polymorphic reference NUMTs
```

`NUMT_detection_results.tsv`

| Column      | Description                                                                                |
| ----------- | ------------------------------------------------------------------------------------------ |
| **numt_id** | Unique identifier for the detected NUMT (e.g., `REFNUMT_1`).                               |
| **source**  | Origin of detection: `reference` (reference genome) or `pangenome` (graph-derived).        |
| **type**    | Classification: `fixed` or `polymorphic`.                                                  |
| **variant** | NUMT-associated variants (for polymorphic NUMTs).                                          |
| **nuchrom** | Nuclear chromosome.                                                                        |
| **nupos**   | Position on the nuclear genome.                                                            |
| **nustart** | Start coordinate on the nuclear genome.                                                    |
| **nuend**   | End coordinate on the nuclear genome.                                                      |
| **mtchrom** | chrM                                                                                       |
| **mtstart** | Start coordinate on the mitochondrial genome.                                              |
| **mtend**   | End coordinate on the mitochondrial genome.                                                |
| **strand**  | Orientation (`+` or `-`).                                                 |

---

## 5. Examples

To test the pipeline with example data:

1. **Prepare data:**

Download the T2T-CHM13v2.0 chromosome 1 FASTA from NCBI and place it in `examples/test_data/`.

2. **Run Snakemake using the example configuration:**

```bash
cd examples
snakemake -p -s ../Snakefile
```

---

## Citation

Coming soon!
