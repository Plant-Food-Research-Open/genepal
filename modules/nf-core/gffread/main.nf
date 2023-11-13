process GFFREAD {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.1--h8b12597_0' :
        'biocontainers/gffread:0.12.1--h8b12597_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.gtf")  , emit: gtf, optional: true
    tuple val(meta), path("*.gff3") , emit: gff, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args         ?: ''
    def prefix      = task.ext.prefix       ?: "${gff.baseName}"
    def extension   = args.contains("-T")   ?  '.gtf' : '.gff3'
    """
    gffread \\
        $gff \\
        $args \\
        -o ${prefix}.${extension}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}
