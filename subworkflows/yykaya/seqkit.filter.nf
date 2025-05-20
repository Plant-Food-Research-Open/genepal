process SEQKIT_GET_LENGTH {
    tag "${meta.id}"
    label 'process_medium'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(genome_fasta)

    output:
    tuple val(meta), path("filtered_${meta.id}.fasta"), path("${meta.id}_contig_list.txt"), emit: filtered_fasta

    script:
    """
    # Filter contigs based on length and output filtered FASTA
    seqkit seq --min-len ${params.min_contig_length} ${genome_fasta} > filtered_${meta.id}.fasta

    # Generate a list of filtered contigs
    seqkit fx2tab --length --name filtered_${meta.id}.fasta | awk '{print \$1}' > ${meta.id}_contig_list.txt
    """
