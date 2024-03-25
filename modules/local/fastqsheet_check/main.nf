// Source:
// https://github.com/nf-core/rnaseq
// MIT: https://github.com/nf-core/rnaseq/blob/master/LICENSE
//
// Changes:
// Added channel permissible_target_assemblies

process FASTQSHEET_CHECK {
    tag "$fastqsheet"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path fastqsheet
    val permissible_target_assemblies

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_fastqsheet.py \\
        $fastqsheet \\
        "$permissible_target_assemblies" \\
        fastqsheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
