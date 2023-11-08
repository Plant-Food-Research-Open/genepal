nextflow.enable.dsl=2

// https://github.com/nf-core/rnaseq
// MIT: https://github.com/nf-core/rnaseq/blob/master/LICENSE
//
// Changes:
// Added channel permissible_target_assemblies
// Changed file name from input_check.nf to extract_samples.nf
// Removed strandedness
// Nowing emitting an extra channel 'assemblies' which indicates the
// assemblies targeted by each read
//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow EXTRACT_SAMPLES {
    take:
    samplesheet                     // file: /path/to/samplesheet.csv
    permissible_target_assemblies   // val: assembly_a,assembly_b

    main:
    SAMPLESHEET_CHECK ( samplesheet, permissible_target_assemblies )
    .csv
    | splitCsv ( header:true, sep:',' )
    | map { create_fastq_channel(it) }
    | set { ch_reads }

    reads = ch_reads.map { meta, fastq -> [[id:meta.id, single_end:meta.single_end], fastq]}
    
    ch_reads
    | flatMap { meta, fastq ->
        meta.target_assemblies.collect { assembly -> [[id:meta.id, single_end:meta.single_end], assembly] }
    }
    | set { assemblies }
    
    emit:
    reads                                       // channel: [ val(meta), [ reads ] ]
    assemblies                                  // channel: [ val(meta), val(assembly) ]
    versions = SAMPLESHEET_CHECK.out.versions   // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id                 = row.sample
    meta.single_end         = row.single_end.toBoolean()
    meta.target_assemblies  = row.target_assemblies.split(";").sort()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}