process BASECALLER {
clusterOptions "--partition cuda --gres=gpu:1"
containerOptions "--nv"
publishDir 'basecall', mode:'copy'

input:
	path pod5_dir
	val kit
output:
	path("dorado_basecalled.bam"), emit:bam
	val 'done', emit:done
script:
"""
dorado basecaller /pixi_env/.pixi/envs/default/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 pod5/ --kit-name ${kit}  > dorado_basecalled.bam
"""
}

process DEMULTIPLEX{
	cpus 8
	publishDir 'demuxed', mode:'copy'
input:
	path bam_file
	path samp_csv
output:
    path("*.fastq.gz"), emit: fastq_gz
	//tuple val(sample_id), path("*.fastq.gz"), emit: fastq_gz
    val 'done', emit: done
script:
"""
mkdir -p excluded/
dorado demux --emit-fastq --output-dir demuxed --no-classify ${bam_file} 
mv demuxed/*/*/*/*/*/*.fastq .

declare -A SAMPLE_MAP
while IFS=';' read -r bc_num sample_name; do
	bc_num=\$(echo "\$bc_num"     | tr -d '[:space:]')
	sample_name=\$(echo "\$sample_name" | tr -d '[:space:]\\r')
	[[ -z "\$bc_num" || -z "\$sample_name" ]] && continue
	SAMPLE_MAP["\$bc_num"]="\$sample_name"
done < ${samp_csv}
for fq in *.fastq; do
	[[ -e "\$fq" ]] || continue          # glob miss guard

	if [[ "\$fq" == *unknown* ]]; then
		mv "\$fq" excluded/
		continue
	fi


	if [[ "\$fq" =~ barcode([0-9]+) ]]; then
		bc_num="\${BASH_REMATCH[1]}"
		bc_num=\$(printf '%02d' "\$((10#\$bc_num))")
	else
		mv "\$fq" excluded/
		continue
	fi

	sample_name="\${SAMPLE_MAP[\$bc_num]:-}"

	if [[ -z "\$sample_name" ]]; then
		mv "\$fq" excluded/
		continue
	fi

	new_name="\${sample_name}-barcode\${bc_num}.fastq"
	mv "\$fq" "\$new_name"
	gzip "\$new_name"          
done
"""
}

process METHYLATION_6mA {
clusterOptions "--partition cuda --gres=gpu:1"
containerOptions "--nv"
publishDir 'methylation', mode:'copy'

input:
	path pod5_dir
output:
	path("*.bam"), emit:bam
	val 'done', emit:done
script:
"""
dorado basecaller /idle/ric.cirillo/zinola.alma/SOFTW/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 pod5/ --modified-bases-models /idle/ric.cirillo/zinola.alma/SOFTW/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1 --kit-name SQK-NBD114-96 > 6mA-met.bam
"""
}

process DEMULTIPLEX_6mA{
	cpus 8
	publishDir 'demuxed_6mA', mode:'copy'
input:
	path bam_file
	path samp_csv
output:
    path("*.bam"), emit: bam
    val 'done', emit: done
script:
"""
mkdir -p excluded
dorado demux --output-dir demuxed_6mA --kit-name SQK-NBD114-96 ${bam_file}
mv demuxed_6mA/*/*/*/*/*/*.bam .

declare -A SAMPLE_MAP
while IFS=';' read -r bc_num sample_name; do
	bc_num=\$(echo "\$bc_num"     | tr -d '[:space:]')
	sample_name=\$(echo "\$sample_name" | tr -d '[:space:]\\r')
	[[ -z "\$bc_num" || -z "\$sample_name" ]] && continue
	SAMPLE_MAP["\$bc_num"]="\$sample_name"
done < ${samp_csv}
for fq in *.bam; do
	[[ -e "\$fq" ]] || continue          

	if [[ "\$fq" == *unknown* ]]; then
		mv "\$fq" excluded/
		continue
	fi

	if [[ "\$fq" =~ barcode([0-9]+) ]]; then
		bc_num="\${BASH_REMATCH[1]}"
		bc_num=\$(printf '%02d' "\$((10#\$bc_num))")
	else
		mv "\$fq" excluded/
		continue
	fi

	sample_name="\${SAMPLE_MAP[\$bc_num]:-}"

	if [[ -z "\$sample_name" ]]; then
		mv "\$fq" excluded/
		continue
	fi

	new_name="\${sample_name}-barcode\${bc_num}.bam"
	mv "\$fq" "\$new_name"
done
"""
}

process MAPPING_6mA{
	conda "/idle/ric.cirillo/common_envs/conda/bcftools"
	cpus 8
	publishDir 'mapped_6mA', mode:'copy'
	tag "$sample_id"
input:
    tuple val(sample_id), path(bam_file)

output:
	tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: bam
	val 'done', emit: done

script:
"""
dorado aligner /beegfs/datasets/buffer/ric.cirillo/MTB/M._tuberculosis_H37Rv_2015-11-13.fasta ${bam_file} | samtools sort -o "${sample_id}.aligned.bam"
samtools index "${sample_id}.aligned.bam"
"""
}

process MODKIT_6mA{
	conda "/idle/ric.cirillo/common_envs/conda/modkit"
	cpus 8
	tag "$sample_id"
	publishDir 'mapped_6mA', mode:'copy'
input:
    tuple val(sample_id), path(bam_file), path(bai_file)

output:
	path("*.bedmethyl"), emit: bedmethyl
	path("*.bedmethyl.log"), emit: log
	val 'done', emit: done

script:
"""
modkit pileup ${bam_file} ${sample_id}.bedmethyl --threads ${task.cpus} --log-filepath ${sample_id}.bedmethyl.log
"""
}

process METHYLATION_4mC {
clusterOptions "--partition cuda --gres=gpu:1"
containerOptions "--nv"
publishDir 'methylation', mode:'copy'

input:
	path pod5_dir
output:
	path("*.bam"), emit:bam
	val 'done', emit:done
script:
"""
dorado basecaller /idle/ric.cirillo/zinola.alma/SOFTW/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 pod5/ --modified-bases-models /idle/ric.cirillo/zinola.alma/SOFTW/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1 > 4mC_5m-met.bam
"""
}

process DEMULTIPLEX_4mC{
	cpus 8
	publishDir 'demuxed_4mC', mode:'copy'
input:
	path bam_file
	path samp_csv
output:
    path("*.bam"), emit: bam
    val 'done', emit: done
script:
"""
mkdir -p excluded
dorado demux --output-dir demuxed_4mC --kit-name SQK-NBD114-96 ${bam_file}
mv demuxed_4mC/*/*/*/*/*/*.bam .

declare -A SAMPLE_MAP
while IFS=';' read -r bc_num sample_name; do
	bc_num=\$(echo "\$bc_num"     | tr -d '[:space:]')
	sample_name=\$(echo "\$sample_name" | tr -d '[:space:]\\r')
	[[ -z "\$bc_num" || -z "\$sample_name" ]] && continue
	SAMPLE_MAP["\$bc_num"]="\$sample_name"
done < ${samp_csv}
for fq in *.bam; do
	[[ -e "\$fq" ]] || continue         

	if [[ "\$fq" == *unknown* ]]; then
		mv "\$fq" excluded/
		continue
	fi

	if [[ "\$fq" =~ barcode([0-9]+) ]]; then
		bc_num="\${BASH_REMATCH[1]}"
		bc_num=\$(printf '%02d' "\$((10#\$bc_num))")
	else
		mv "\$fq" excluded/
		continue
	fi

	sample_name="\${SAMPLE_MAP[\$bc_num]:-}"

	if [[ -z "\$sample_name" ]]; then
		mv "\$fq" excluded/
		continue
	fi

	new_name="\${sample_name}-barcode\${bc_num}.bam"
	mv "\$fq" "\$new_name"
done
"""
}

process MAPPING_4mC{
	conda "/idle/ric.cirillo/common_envs/conda/bcftools"
	cpus 8
	publishDir 'mapped_4mC', mode:'copy'
	tag "$sample_id"
input:
    tuple val(sample_id), path(bam_file)

output:
	tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: bam
	val 'done', emit: done

script:
"""
dorado aligner /beegfs/datasets/buffer/ric.cirillo/MTB/M._tuberculosis_H37Rv_2015-11-13.fasta ${bam_file} | samtools sort -o "${sample_id}.aligned.bam"
samtools index "${sample_id}.aligned.bam"
"""
}

process MODKIT_4mC{
	conda "/idle/ric.cirillo/common_envs/conda/modkit"
	cpus 8
	tag "$sample_id"
	publishDir 'mapped_4mC', mode:'copy'
input:
    tuple val(sample_id), path(bam_file), path(bai_file)

output:
	path("*.bedmethyl"), emit: bedmethyl
	path("*.bedmethyl.log"), emit: log
	val 'done', emit: done

script:
"""
modkit pileup ${bam_file} ${sample_id}.bedmethyl --threads ${task.cpus} --log-filepath ${sample_id}.bedmethyl.log
"""
}

process METHYLATION_5mCG {
clusterOptions "--partition cuda --gres=gpu:1"
containerOptions "--nv"
publishDir 'methylation', mode:'copy'

input:
	path pod5_dir
output:
	path("*.bam"), emit:bam
	val 'done', emit:done
script:
"""
dorado basecaller /idle/ric.cirillo/zinola.alma/SOFTW/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 pod5/ --modified-bases-models /idle/ric.cirillo/zinola.alma/SOFTW/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2 > 5mCG_5hm-met.bam
"""
}

process DEMULTIPLEX_5mCG{
	cpus 8
	publishDir 'demuxed_5mCG', mode:'copy'
input:
	path bam_file
	path samp_csv
output:
    path("*.bam"), emit: bam
    val 'done', emit: done
script:
"""
mkdir -p excluded
dorado demux --output-dir demuxed_5mCG --kit-name SQK-NBD114-96 ${bam_file}
mv demuxed_5mCG/*/*/*/*/*/*.bam .

declare -A SAMPLE_MAP

while IFS=';' read -r bc_num sample_name; do
	bc_num=\$(echo "\$bc_num"     | tr -d '[:space:]')
	sample_name=\$(echo "\$sample_name" | tr -d '[:space:]\\r')
	[[ -z "\$bc_num" || -z "\$sample_name" ]] && continue
	SAMPLE_MAP["\$bc_num"]="\$sample_name"
done < ${samp_csv}
for fq in *.bam; do
	[[ -e "\$fq" ]] || continue         

	if [[ "\$fq" == *unknown* ]]; then
		mv "\$fq" excluded/
		continue
	fi

	if [[ "\$fq" =~ barcode([0-9]+) ]]; then
		bc_num="\${BASH_REMATCH[1]}"
		bc_num=\$(printf '%02d' "\$((10#\$bc_num))")
	else
		mv "\$fq" excluded/
		continue
	fi

	sample_name="\${SAMPLE_MAP[\$bc_num]:-}"

	if [[ -z "\$sample_name" ]]; then
		mv "\$fq" excluded/
		continue
	fi

	new_name="\${sample_name}-barcode\${bc_num}.bam"
	mv "\$fq" "\$new_name"
done
"""
}

process MAPPING_5mCG{
	conda "/idle/ric.cirillo/common_envs/conda/bcftools"
	cpus 8
	publishDir 'mapped_5mCG', mode:'copy'
	tag "$sample_id"
input:
    tuple val(sample_id), path(bam_file)

output:
	tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: bam
	val 'done', emit: done

script:
"""
dorado aligner /beegfs/datasets/buffer/ric.cirillo/MTB/M._tuberculosis_H37Rv_2015-11-13.fasta ${bam_file} | samtools sort -o "${sample_id}.aligned.bam"
samtools index "${sample_id}.aligned.bam"
"""
}

process MODKIT_5mCG{
	conda "/idle/ric.cirillo/common_envs/conda/modkit"
	cpus 8
	tag "$sample_id"
	publishDir 'mapped_5mCG', mode:'copy'
input:
    tuple val(sample_id), path(bam_file), path(bai_file)

output:
	path("*.bedmethyl"), emit: bedmethyl
	path("*.bedmethyl.log"), emit: log
	val 'done', emit: done

script:
"""
modkit pileup ${bam_file} ${sample_id}.bedmethyl --threads ${task.cpus} --log-filepath ${sample_id}.bedmethyl.log
"""
}

process METHYLATION_5mC {
clusterOptions "--partition cuda --gres=gpu:1"
containerOptions "--nv"
publishDir 'methylation', mode:'copy'

input:
	path pod5_dir
output:
	path("*.bam"), emit:bam
	val 'done', emit:done
script:
"""
dorado basecaller /idle/ric.cirillo/zinola.alma/SOFTW/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 pod5/ --modified-bases-models /idle/ric.cirillo/zinola.alma/SOFTW/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mC_5hmC@v2 > 5mC_5hmC-met.bam
"""
}

process DEMULTIPLEX_5mC{
	cpus 8
	publishDir 'demuxed_5mC', mode:'copy'
input:
	path bam_file
	path samp_csv
output:
    path("*.bam"), emit: bam
    val 'done', emit: done
script:
"""
mkdir -p excluded
dorado demux --output-dir demuxed_5mC --kit-name SQK-NBD114-96 ${bam_file}
mv demuxed_5mC/*/*/*/*/*/*.bam .

declare -A SAMPLE_MAP
while IFS=';' read -r bc_num sample_name; do
	bc_num=\$(echo "\$bc_num"     | tr -d '[:space:]')
	sample_name=\$(echo "\$sample_name" | tr -d '[:space:]\\r')
	[[ -z "\$bc_num" || -z "\$sample_name" ]] && continue
	SAMPLE_MAP["\$bc_num"]="\$sample_name"
done < ${samp_csv}
for fq in *.bam; do
	[[ -e "\$fq" ]] || continue          

	# exclude unknown sample
	if [[ "\$fq" == *unknown* ]]; then
		mv "\$fq" excluded/
		continue
	fi

	# extract barcode
	# Matches "barcode" followed by digits, e.g. barcode01, barcode12
	if [[ "\$fq" =~ barcode([0-9]+) ]]; then
		bc_num="\${BASH_REMATCH[1]}"
		bc_num=\$(printf '%02d' "\$((10#\$bc_num))")
	else
		mv "\$fq" excluded/
		continue
	fi

	# sample name
	sample_name="\${SAMPLE_MAP[\$bc_num]:-}"

	if [[ -z "\$sample_name" ]]; then
		# barcode exists but has no entry in CSV -> exclude
		mv "\$fq" excluded/
		continue
	fi

	# rename + zip
	new_name="\${sample_name}-barcode\${bc_num}.bam"
	mv "\$fq" "\$new_name"
done
"""
}

process MAPPING_5mC{
	conda "/idle/ric.cirillo/common_envs/conda/bcftools"
	cpus 8
	publishDir 'mapped_5mC', mode:'copy'
	tag "$sample_id"
input:
    tuple val(sample_id), path(bam_file)

output:
	tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: bam
	val 'done', emit: done

script:
"""
dorado aligner /beegfs/datasets/buffer/ric.cirillo/MTB/M._tuberculosis_H37Rv_2015-11-13.fasta ${bam_file} | samtools sort -o "${sample_id}.aligned.bam"
samtools index "${sample_id}.aligned.bam"
"""
}

process MODKIT_5mC{
	conda "/idle/ric.cirillo/common_envs/conda/modkit"
	cpus 8
	tag "$sample_id"
	publishDir 'mapped_5mC', mode:'copy'
input:
    tuple val(sample_id), path(bam_file), path(bai_file)

output:
	path("*.bedmethyl"), emit: bedmethyl
	path("*.bedmethyl.log"), emit: log
	val 'done', emit: done

script:
"""
modkit pileup ${bam_file} ${sample_id}.bedmethyl --threads ${task.cpus} --log-filepath ${sample_id}.bedmethyl.log
"""
}

process NANOCOMP{
	cpus 10
	publishDir 'Nanocomp', mode:'copy'
input:
	path fastq_gz
output:
	//path("demuxed"), emit: demux_dir
    path("*.txt"), emit: NANOCOMP
    val 'done', emit: done
script:
"""
NanoComp --threads ${task.cpus} --tsv_stats --fastq ${fastq_gz}
"""
}

process KRAKEN{
	cpus 16
	memory '150GB'
	tag "$replicateId"
	publishDir 'kraken', mode:'copy'
input:
    tuple val(replicateId), path(fastq_gz)
	val kraken_db

output:
    tuple val(replicateId), path("${replicateId}.kreport"), emit: kraken
    val 'done', emit: done
script:
"""
mkdir kraken
kraken2 --db /beegfs/datasets/buffer/ric.cirillo/kraken_db/${kraken_db} --threads ${task.cpus} --use-names --gzip-compressed --output kraken/${replicateId}.kraken --report kraken/${replicateId}.kreport ${fastq_gz}
mv kraken/*.kreport .
"""
}

process BRACKEN {
cpus 16
tag "$replicateId"
publishDir "bracken", mode:"copy"

input:
    tuple val(replicateId), path(kreport)
	val kraken_db
output:
    tuple val(replicateId), path("*.bout"), emit : bout
	tuple val(replicateId),path("*.report"), emit: breport
	val 'done', emit:done
script:
"""
mkdir bracken
bracken -d /beegfs/datasets/buffer/ric.cirillo/kraken_db/${kraken_db} -i $kreport -o bracken/${replicateId}.bout -w bracken/${replicateId}.report -r 150
mv bracken/* .
"""
}

process COVERAGE {
    cpus 2
    publishDir 'coverage', mode: 'copy'

input:
	path fastq_gz
	val genome_size

output:
	path "coverage.tsv", emit: coverage_table

script:
"""
GENOME_SIZE=${genome_size}

echo -e "sample\tcoverage" > coverage.tsv

for f in ${fastq_gz}; do

	#sample=\$(basename \$f .fastq.gz)

	bases=\$(zcat \$f | awk 'NR%4==2 {total += length(\$0)} END {print total}')

	coverage=\$(echo "scale=4; \$bases / \$GENOME_SIZE" | bc)

	echo -e "\$f\t\$coverage" >> coverage.tsv
done
"""
}


process FINAL_REPORT {

    publishDir 'output', mode: 'copy'

    input:
        path nanostats
        path coverage_table
        path bracken_bout

    output:
        path "Final_Report.tsv"

    script:
    """
    set -euo pipefail

    NANOSTATS=${nanostats}
    COVERAGE=${coverage_table}
    BRACKEN_FILES="${bracken_bout}"
    OUT=Final_Report.tsv

    # header nanocomp
    read header < \$NANOSTATS
    echo "\$header" > "\$OUT"

    # nanocomp stats
    tail -n +2 \$NANOSTATS >> "\$OUT"

    # sample coverage
    awk '
    NR==FNR { cov[\$1]=\$2; next }
    FNR==1 {
        printf "expected_coverage"
        for (i=2; i<=NF; i++) {
            printf "\\t%s", cov[\$i]
        }
        printf "\\n"
    }' \$COVERAGE <(echo "\$header") >> "\$OUT"

    # Link samples between nanocomp and bracken
    declare -A BRACKEN_MAP
    for f in \$BRACKEN_FILES; do
        bc=\$(basename "\$f" | sed -E 's/.*(barcode[0-9]+).*/\\1/')
        BRACKEN_MAP["\$bc"]="\$f"
    done

    # Major specie bracken
    printf "Major_Species" >> "\$OUT"

    for sample in \$(echo "\$header" | tr '\\t' '\\n' | tail -n +2); do
        bc=\$(echo "\$sample" | sed -E 's/.*(barcode[0-9]+).*/\\1/')
        f="\${BRACKEN_MAP[\$bc]}"

        if [[ -f "\$f" ]]; then
            sp=\$(awk 'BEGIN{FS="\\t"} { if (\$7+0 > max) { max=\$7; name=\$1 } } END{print name}' "\$f")
        else
            sp="NA"
        fi

        printf "\t%s" "\$sp" >> "\$OUT"
    done
    printf "\n" >> "\$OUT"

    # Percentage specie
    printf "Percentage" >> "\$OUT"

    for sample in \$(echo "\$header" | tr '\\t' '\\n' | tail -n +2); do
        bc=\$(echo "\$sample" | sed -E 's/.*(barcode[0-9]+).*/\\1/')
        f="\${BRACKEN_MAP[\$bc]}"

        if [[ -f "\$f" ]]; then
            pct=\$(awk 'BEGIN{FS="\\t"} {if (\$7+0 > max){max=\$7}} END{print max}' "\$f")
        else
            pct="NA"
        fi

        printf "\t%s" "\$pct" >> "\$OUT"
    done
    printf "\n" >> "\$OUT"
    """
}
