# nf_ont:
This pipeline is designed to process Oxford Nanopore Technologies (ONT) sequencing reads and generate a summary report for selected samples. The report includes expected sequencing coverage and Bracken-based taxonomic analysis.

# Input Requirement:

The pipeline requires the following inputs:
- A directory named pod5/ containing all .pod5 raw signal files.
- A sample sheet provided via the --file parameter.
The sample sheet must list samples and their associated barcodes using a semicolon (;) as a separator.
For example:
01;SAMPA
02;SAMPB
The files .pod5 that are not associated to a barcode will be put in the directory excluded/

# Optional analysis:
The pipeline supports optional DNA methylation analysis. The following models are available:
- 6mA
- 4mC
- 5mC
- 5mCG
You can choose which methylation analysis to run. By default, this option is set to false.

The pipeline can optionally perform:
- Read alignment using the Dorado aligner
- Methylation extraction using Modkit pileup
By default, this step is disabled (false).
