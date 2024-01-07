// Source:
// https://github.com/nf-core/rnaseq
// MIT: https://github.com/nf-core/rnaseq/blob/master/LICENSE
//
// Check input samplesheet and get read channels
//
// Changes:
// Added channel permissible_target_assemblies
// Changed file name from input_check.nf to extract_samples.nf
// Removed strandedness
// Nowing emitting an extra channel 'assemblies' which indicates the
// assemblies targeted by each read

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow EXTRACT_SAMPLES {
    take:
    samplesheet                     // file: /path/to/samplesheet.csv
    permissible_target_assemblies   // val: assembly_a,assembly_b

    main:
    SAMPLESHEET_CHECK ( samplesheet, permissible_target_assemblies )
    .csv
    | splitCsv ( header:true, sep:',' )
    | combine ( samplesheet )
    | map { row, sheet ->
        create_fastq_channel(row, sheet.getParent())
    }
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
def create_fastq_channel(LinkedHashMap row, sheetPath) {
    // create meta map
    def meta = [:]
    meta.id                 = row.sample
    meta.single_end         = row.single_end.toBoolean()
    meta.target_assemblies  = row.target_assemblies.split(";").sort()

    def fq1                 = row.fastq_1.startsWith("/") ? row.fastq_1 : "$sheetPath/${row.fastq_1}"
    def fq2                 = row.fastq_2.startsWith("/") ? row.fastq_2 : "$sheetPath/${row.fastq_2}"

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(fq1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${fq1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(fq1) ] ]
    } else {
        if (!file(fq2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${fq2}"
        }
        fastq_meta = [ meta, [ file(fq1), file(fq2) ] ]
    }

    println fastq_meta

    return fastq_meta
}