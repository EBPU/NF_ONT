# NF_ONT

## Overview

**NF_ONT** is a Nextflow pipeline designed for the analysis of Oxford Nanopore Technologies (ONT) sequencing data. The pipeline processes raw POD5 files, performs basecalling and demultiplexing, and generates a report for selected samples.

The standard output includes:

* Sequencing coverage estimation
* Taxonomic profiling using Bracken
* Quality control metrics

Optionally, the pipeline can perform DNA methylation analysis using Dorado and Modkit, generating dedicated methylation reports for the selected modification types.

---

## Input Requirements

The pipeline requires the following inputs:

### 1. POD5 Directory

A directory containing all raw ONT signal files in `.pod5` format:

```text
pod5/
├── file1.pod5
├── file2.pod5
└── ...
```

### 2. Sample Sheet

A sample sheet provided through the `params.list` parameter.

The file must contain barcode-to-sample associations separated by semicolons (`;`):

```text
01;SAMPA
02;SAMPB
03;SAMPC
```

Reads assigned to barcodes that are not listed in the sample sheet will be automatically placed in the `excluded/` directory.

---

## Optional DNA Methylation Analysis

The pipeline supports optional methylation calling using Dorado models.

Available methylation models:

* `6mA`
* `4mC`
* `5mC`
* `5mCG`

Each methylation type can be enabled independently. By default, all methylation analyses are disabled.

### Optional Alignment and Methylation Extraction

For each selected methylation model, the pipeline can additionally perform:

* Read alignment using the Dorado aligner
* Methylation extraction and aggregation using Modkit (`pileup`)

These steps are disabled by default.

---

## Reference Genome

The default reference genome is:

**Mycobacterium tuberculosis H37Rv**

To use a custom reference:

1. Place the FASTA file inside the `REF/` directory.
2. Specify the reference using the `--ref` parameter.
3. Update the genome size using `--genome_size`.

---

## Parameters

| Parameter      | Description                                         | Default                            |
| -------------- | --------------------------------------------------- | ---------------------------------- |
| `params.pod5`         | Directory containing POD5 files                     | `pod5/`                            |
| `params.kit`          | ONT barcoding kit used for demultiplexing           | `SQK-NBD114-96`                    |
| `params.list`         | Sample sheet containing barcode/sample associations | `samples.csv`                      |
| `params.ref`          | Reference genome used for analysis                  | `Mycobacterium tuberculosis H37Rv` |
| `params.run_6mA`      | Enable 6mA methylation calling                      | `false`                            |
| `params.run_4mC`      | Enable 4mC methylation calling                      | `false`                            |
| `params.run_5mCG`     | Enable 5mCG methylation calling                     | `false`                            |
| `params.run_5mC`      | Enable 5mC methylation calling                      | `false`                            |
| `params.mapping_6mA`  | Run alignment and Modkit analysis for 6mA calls     | `false`                            |
| `params.mapping_4mC`  | Run alignment and Modkit analysis for 4mC calls     | `false`                            |
| `params.mapping_5mCG` | Run alignment and Modkit analysis for 5mCG calls    | `false`                            |
| `params.mapping_5mC`  | Run alignment and Modkit analysis for 5mC calls     | `false`                            |

---

## Usage

### Standard Analysis

Run the pipeline using the default reference genome:

```bash
nextflow run NF_ONT
```

### 6mA Methylation Analysis

Run 6mA methylation calling together with alignment and Modkit processing:

```bash
nextflow run NF_ONT --run_6mA true --mapping_6mA true
```

### Custom Reference Genome

```bash
nextflow run NF_ONT --ref custom_reference.fasta --genome_size 4500000
```

---

## Output

The pipeline generates:

* Demultiplexed FASTQ files
* Coverage statistics
* Bracken taxonomic classification reports
* Quality control summaries
* Optional methylation reports (Dorado + Modkit)
* Final aggregated sample reports

