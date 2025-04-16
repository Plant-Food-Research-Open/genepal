process CUSTOM_RMOUTTOGFF3 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-bioperl:1.7.8--hdfd78af_1':
        'biocontainers/perl-bioperl:1.7.8--hdfd78af_1' }"

    input:
    tuple val(meta), path(rmout)

    output:
    tuple val(meta), path("*.gff3") , emit: gff3
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'rmouttogff3.pl'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | sed -n 's|This is perl.*(v\\(.*\\)) .*|\\1|p'  )
    END_VERSIONS
    """
}
