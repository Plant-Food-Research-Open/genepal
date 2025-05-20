process SEQKIT {
    tag "$sample_id"
    label 'seqkit'
    
    input:
    tuple val(sample_id), path(contigs)
    val min_length
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.fasta")
    
    script:
    """
    seqkit seq --min-len ${min_length} ${contigs} > ${sample_id}.filtered.fasta
    """
}
