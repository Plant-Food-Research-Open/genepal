process BRAKER3 {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::braker3=3.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'registry.hub.docker.com/teambraker/braker3:v.1.0.4':
        'registry.hub.docker.com/teambraker/braker3:v.1.0.4' }"

    input:
    tuple val(meta), path(fasta), path(rnaseq_sets_ids), path(rnaseq_sets_dirs), path(bam), path(proteins), path(hintsfile)

    output:
    tuple val(meta), path("${prefix}/braker.gtf")      , emit: gtf
    tuple val(meta), path("${prefix}/braker.codingseq"), emit: cds
    tuple val(meta), path("${prefix}/braker.aa")       , emit: aa
    tuple val(meta), path("${prefix}/hintsfile.gff")   , emit: hintsfile
    tuple val(meta), path("${prefix}/braker.log")      , emit: log
    tuple val(meta), path("${prefix}/what-to-cite.txt"), emit: citations
    tuple val(meta), path("${prefix}/braker.gff3")     , emit: gff3     , optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def rna_ids  = rnaseq_sets_ids ? "--rnaseq_sets_ids=${rnaseq_sets_ids}" : ''
    def rna_dirs = rnaseq_sets_dirs ? "--rnaseq_sets_dirs=${rnaseq_sets_dirs}" : ''
    def bam      = bam ? "--bam=${bam}" : ''
    def proteins = proteins ? "--prot_seq=${proteins}" : ''
    def hints    = hintsfile ? "--hints=${hintsfile}" : ''
    """
    braker.pl \\
        --genome ${fasta} \\
        --species ${prefix} \\
        --workingdir ${prefix} \\
        --threads ${task.cpus} \\
        ${rna_ids} \\
        ${rna_dirs} \\
        ${bam} \\
        ${proteins} \\
        ${hints} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         braker3: \$(braker.pl --version 2>&1 | sed 's/^.*BRAKER3 v//; s/ .*\$//')
    END_VERSIONS
    """
}
