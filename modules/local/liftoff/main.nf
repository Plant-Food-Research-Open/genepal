process LIFTOFF {
    tag "$meta.id"
    label "process_high"

    container "https://depot.galaxyproject.org/singularity/liftoff:1.6.3--pyhdfd78af_0"

    input:
    tuple val(meta), path(target_fa)
    path ref_fa
    path ref_gff
    
    output:
    tuple val(meta), path("*.liftoff.gff3")         , emit: gff3
    tuple val(meta), path("unmapped_features.txt")  , emit: unmapped
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    liftoff \\
    -g $ref_gff \\
    -p $task.cpus \\
    $args \\
    $target_fa \\
    $ref_fa \\
    > "${prefix}.liftoff.gff3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftoff: \$(liftoff --version)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.liftoff.gff3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftoff: \$(liftoff --version)
    END_VERSIONS
    """
}