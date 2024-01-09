process BRAKER3 {
    tag "${meta.id}"
    label 'process_high'

    container "gallvp/teambraker_braker3:v1.0.6"

    input:
    tuple val(meta), path(fasta)
    path bam
    path rnaseq_sets_dirs
    path rnaseq_sets_ids
    path proteins
    path hintsfile

    output:
    tuple val(meta), path("${prefix}/braker.gtf")      , emit: gtf
    tuple val(meta), path("${prefix}/braker.codingseq"), emit: cds
    tuple val(meta), path("${prefix}/braker.aa")       , emit: aa
    tuple val(meta), path("${prefix}/hintsfile.gff")   , emit: hintsfile, optional: true
    tuple val(meta), path("${prefix}/braker.log")      , emit: log
    tuple val(meta), path("${prefix}/what-to-cite.txt"), emit: citations
    tuple val(meta), path("${prefix}/braker.gff3")     , emit: gff3     , optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args                                     ?: ''
    prefix          = task.ext.prefix                                   ?: "${meta.id}"

    def test_mode   = args.contains('--testMode') // Custom flag for test data
    def args_fmt    = test_mode ? args.replace('--testMode', '') : args

    def rna_ids     = rnaseq_sets_ids           ? "--rnaseq_sets_ids=${rnaseq_sets_ids}"    : ''
    def rna_dirs    = rnaseq_sets_dirs          ? "--rnaseq_sets_dirs=${rnaseq_sets_dirs}"  : ''
    def bam         = bam && !test_mode         ? "--bam=${bam}"                            : ''
    def proteins    = proteins && !test_mode    ? "--prot_seq=${proteins}"                  : ''
    def hints       = hintsfile                 ? "--hints=${hintsfile}"                    : ''

    def new_species = args.contains('--species')   ? '' : "--species new_species"
    """
    cp -r /usr/share/augustus/config augustus_config

    braker.pl \\
        --genome ${fasta} \\
        ${new_species} \\
        --workingdir ${prefix} \\
        --AUGUSTUS_CONFIG_PATH "\$(pwd)/augustus_config" \\
        --threads ${task.cpus} \\
        ${rna_ids} \\
        ${rna_dirs} \\
        ${bam} \\
        ${proteins} \\
        ${hints} \\
        ${args_fmt}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker3: \$(braker.pl --version 2>/dev/null | sed 's/braker.pl version//')
    END_VERSIONS
    """

    stub:
    prefix          = task.ext.prefix                       ?: "${meta.id}"
    def rna_ids     = rnaseq_sets_ids                       ? "--rnaseq_sets_ids=${rnaseq_sets_ids}"    : ''
    def touch_hints = (rna_ids || bam || proteins || hints) ? "touch ${prefix}/hintsfile.gff"           : ''
    """
    mkdir "$prefix"

    touch "${prefix}/braker.gtf"
    touch "${prefix}/braker.codingseq"
    touch "${prefix}/braker.aa"
    $touch_hints
    touch "${prefix}/braker.log"
    touch "${prefix}/what-to-cite.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker3: \$(braker.pl --version 2>/dev/null | sed 's/braker.pl version//')
    END_VERSIONS
    """
}
