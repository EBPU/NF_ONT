nextflow.enable.dsl = 2
/*
 * Define the default parameters
 */ 
	params.pod5	= "$baseDir/pod5/"
	params.results	= "OUTPUT"
	params.ref = "$baseDir"
	params.list = "$baseDir/samples.csv"
    params.run_6mA = false
    params.run_4mC = false
    params.run_5mCG = false
    params.run_5mC = false
    params.genome_size = 4411529
    kraken_db = "standard_db"
    params.mapping_6mA = false
    params.mapping_4mC = false
    params.mapping_5mCG = false
    params.mapping_5mC = false

include {
    BASECALLER;
    DEMULTIPLEX;
    METHYLATION_6mA;
    DEMULTIPLEX_6mA;
    MAPPING_6mA;
    MODKIT_6mA;
    METHYLATION_4mC;
    DEMULTIPLEX_4mC;
    MAPPING_4mC;
    MODKIT_4mC;
    METHYLATION_5mCG;
    DEMULTIPLEX_5mCG;
    MAPPING_5mCG;
    MODKIT_5mCG;
    METHYLATION_5mC;
    DEMULTIPLEX_5mC;
    MAPPING_5mC;
    MODKIT_5mC;
    NANOCOMP;
    COVERAGE;
    KRAKEN;
    BRACKEN;
    FINAL_REPORT
	} from "$baseDir/module.nf"
workflow {
    samp_ch = Channel.fromPath(params.list)
    basecaller_out = BASECALLER(params.pod5)
    demultiplex_out = DEMULTIPLEX(basecaller_out.bam,samp_ch)
    reads_with_id = demultiplex_out.fastq_gz
    .flatten()
    .map { file ->
        def sample = file.getName().replaceAll(/\.fastq\.gz$/, '')
        tuple(sample, file)
    }
    nanocomp_out = NANOCOMP(demultiplex_out.fastq_gz)
    KRAKEN(reads_with_id,params.kraken_db)
    bracken_out = BRACKEN(KRAKEN.out.kraken,params.kraken_db)
    coverage_out = COVERAGE(demultiplex_out.fastq_gz.collect(), params.genome_size)
    FINAL_REPORT(nanocomp_out.NANOCOMP,coverage_out,bracken_out.bout.map { it[1] }.collect())
    if (params.run_6mA){
        basecaller_out_6mA=METHYLATION_6mA(params.pod5)
        demux_6mA_out = DEMULTIPLEX_6mA(basecaller_out_6mA.bam,samp_ch) 
        if (params.mapping_6mA){
            mapped_6mA_out = MAPPING_6mA(demux_6mA_out.bam.flatten().map { file -> tuple(file.baseName, file)})
            MODKIT_6mA(mapped_6mA_out.bam)
        }
    }
    if (params.run_4mC){
        basecaller_out_4mC=METHYLATION_4mC(params.pod5)
        DEMULTIPLEX_4mC(basecaller_out_4mC.bam,samp_ch)
        if (params.mapping_4mC){
            mapped_4mC_out = MAPPING_4mC(demux_4mC_out.bam.flatten().map { file -> tuple(file.baseName, file)})
            MODKIT_4mC(mapped_4mC_out.bam)
        }
    }
    if (params.run_5mCG){
        basecaller_out_5mCG=METHYLATION_5mCG(params.pod5)
        DEMULTIPLEX_5mCG(basecaller_out_5mCG.bam,samp_ch)
        if (params.mapping_5mCG){
            mapped_5mCG_out = MAPPING_5mCG(demux_5mCG_out.bam.flatten().map { file -> tuple(file.baseName, file)})
            MODKIT_5mCG(mapped_5mCG_out.bam)
        }
    }
    if (params.run_5mC){
        basecaller_out_5mC=METHYLATION_5mC(params.pod5)
        DEMULTIPLEX_5mC(basecaller_out_5mC.bam,samp_ch)
        if (params.mapping_5mC){
            mapped_5mC_out = MAPPING_5mC(demux_5mC_out.bam.flatten().map { file -> tuple(file.baseName, file)})
            MODKIT_5mC(mapped_5mC_out.bam)
        }
    }
}