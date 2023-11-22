process REPEATMASKER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::repeatmasker=4.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(lib)

    output:
    tuple val(meta), path("${meta.id}/*.f*a.masked") , emit: fasta_masked
    tuple val(meta), path("${meta.id}/*.f*a.out")    , emit: fasta_out
    tuple val(meta), path("${meta.id}/*.f*a.tbl")    , emit: fasta_tbl
    tuple val(meta), path("${meta.id}/*.f*a.cat.gz") , emit: fasta_cat_gz , optional: true
    tuple val(meta), path("${meta.id}/*.f*a.out.gff"), emit: fasta_out_gff, optional: true
    tuple val(meta), path("${meta.id}/*.f*a.align")  , emit: fasta_align  , optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.1.5'  // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    RepeatMasker \\
        -lib ${lib} \\
        -pa ${task.cpus} \\
        -dir ${prefix} \\
        ${args} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: ${VERSION}
    END_VERSIONS
    """
}
