process FASTAVALIDATE {
    tag "$meta.id"
    label 'process_single'

    // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker://gallvp/fasta_validator:a6a2ec1_ps'

    input:
    tuple val(meta), path(fasta)
    
    output:
    tuple val(meta), path('*.validated.fasta')  , emit: valid_fasta , optional: true
    tuple val(meta), path('*.error.log')        , emit: error_log   , optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fasta_validate \\
        -v $fasta \\
        2> "${prefix}.error.log" \\
        || echo "Errors from fasta_validate printed to ${prefix}.error.log"
    
    if [ \$(cat "${prefix}.error.log" | wc -l) -gt 0 ]; then
        echo "Validation failed..."
        cat "${prefix}.error.log"
    else
        rm "${prefix}.error.log"
        
        cat $fasta \\
            > "${prefix}.validated.fasta"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasta_validate: \$(md5sum \$(which fasta_validate) | cut -d' ' -f1)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.validated.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasta_validate: \$(md5sum \$(which fasta_validate) | cut -d' ' -f1)
    END_VERSIONS
    """
}
