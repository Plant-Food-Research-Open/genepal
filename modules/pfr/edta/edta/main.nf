process EDTA_EDTA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/edta:2.1.0--hdfd78af_1':
        'biocontainers/edta:2.1.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)
    path cds
    path curatedlib
    path rmout
    path exclude

    output:
    tuple val(meta), path('*.log')              , emit: log
    tuple val(meta), path('*.EDTA.TElib.fa')    , emit: te_lib_fasta
    tuple val(meta), path('*.EDTA.pass.list')   , emit: pass_list           , optional: true
    tuple val(meta), path('*.EDTA.out')         , emit: out_file            , optional: true
    tuple val(meta), path('*.EDTA.TEanno.gff3') , emit: te_anno_gff3        , optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args         ?: ''
    def prefix          = task.ext.prefix       ?: "${meta.id}"
    def mod_file_name   = "${fasta}.mod"
    def cds_file        = cds                   ? "--cds $cds"              : ''
    def curatedlib_file = curatedlib            ? "--curatedlib $curatedlib": ''
    def rmout_file      = rmout                 ? "--rmout $rmout"          : ''
    def exclude_file    = exclude               ? "--exclude $exclude"      : ''
    """
    EDTA.pl \\
        --genome $fasta \\
        --threads $task.cpus \\
        $cds_file \\
        $curatedlib_file \\
        $rmout_file \\
        $exclude_file \\
        $args \\
        &> >(tee "${prefix}.log" 2>&1)

    mv \\
        "${mod_file_name}.EDTA.TElib.fa" \\
        "${prefix}.EDTA.TElib.fa"

    [ -f "${mod_file_name}.EDTA.raw/LTR/${mod_file_name}.pass.list" ] \\
        && mv \\
            "${mod_file_name}.EDTA.raw/LTR/${mod_file_name}.pass.list" \\
            "${prefix}.EDTA.pass.list" \\
        || echo "EDTA did not produce a pass.list file"

    [ -f "${mod_file_name}.EDTA.anno/${mod_file_name}.out" ] \\
        && mv \\
            "${mod_file_name}.EDTA.anno/${mod_file_name}.out" \\
            "${prefix}.EDTA.out" \\
        || echo "EDTA did not produce an out file"

    [ -f "${mod_file_name}.EDTA.TEanno.gff3" ] \\
        && mv \\
            "${mod_file_name}.EDTA.TEanno.gff3" \\
            "${prefix}.EDTA.TEanno.gff3" \\
        || echo "EDTA did not produce a TEanno gff3 file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EDTA: \$(EDTA.pl -h | awk ' /##### Extensive/ {print \$7}')
    END_VERSIONS
    """

    stub:
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def touch_pass_list = args.contains("--anno 1") ? "touch ${prefix}.EDTA.pass.list"  : ''
    def touch_out_file  = args.contains("--anno 1") ? "touch ${prefix}.EDTA.out"        : ''
    def touch_te_anno   = args.contains("--anno 1") ? "touch ${prefix}.EDTA.TEanno.gff3": ''
    """
    touch "${prefix}.log"
    touch "${prefix}.EDTA.TElib.fa"
    $touch_pass_list
    $touch_out_file
    $touch_te_anno

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EDTA: \$(EDTA.pl -h | awk ' /##### Extensive/ {print \$7}')
    END_VERSIONS
    """
}
