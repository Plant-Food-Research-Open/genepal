process AGAT_SPCOMPLEMENTANNOTATIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0'
        : 'biocontainers/agat:1.4.2--pl5321hdfd78af_0'}"

    input:
    tuple val(meta), path(ref_gff)
    tuple val(meta2), path(add_gffs)
    path config

    output:
    tuple val(meta), path("*.gff"), emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def config_param = config ? "--config ${config}" : ''
    def file_names = "${add_gffs}".split(' ')
    def add_gff_param = file_names.collect { "--add ${it}" }.join(' ')
    if (file_names.contains("${prefix}.gff") || "${ref_gff}" == "${prefix}.gff") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    agat_sp_complement_annotations.pl \\
        --ref ${ref_gff} \\
        ${add_gff_param} \\
        ${config_param} \\
        ${args} \\
        --output ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_complement_annotations.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_names = "${add_gffs}".split(' ')
    if (file_names.contains("${prefix}.gff") || "${ref_gff}" == "${prefix}.gff") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    touch ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_complement_annotations.pl -h | sed -n 's/.*(AGAT) - Version: \\(.*\\) .*/\\1/p')
    END_VERSIONS
    """
}
