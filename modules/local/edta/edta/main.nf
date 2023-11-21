process EDTA {
    tag "$meta.id"
    label "process_high"
    label "process_week_long"
    
    container 'https://depot.galaxyproject.org/singularity/edta:2.1.0--hdfd78af_1'

    input:
    tuple val(meta), path(fasta_file)
    
    output:
    tuple val(meta), path('*.EDTA.TElib.fa')    , emit: te_lib_fasta
    tuple val(meta), path('*.EDTA.intact.gff3') , emit: intact_gff3
    tuple val(meta), path('*.EDTA.pass.list')   , emit: pass_list
    tuple val(meta), path('*.EDTA.out')         , emit: out_file
    tuple val(meta), path('*.EDTA.TEanno.gff3') , emit: te_anno_gff3
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def modFileName = "${fasta_file}.mod"
    """
    EDTA.pl \\
    --genome $fasta_file \\
    --threads $task.cpus \\
    $args
    
    if [ -f "${modFileName}.EDTA.raw/LTR/${modFileName}.pass.list" ]; then
        cat "${modFileName}.EDTA.raw/LTR/${modFileName}.pass.list" \\
        > "${modFileName}.EDTA.pass.list"
    else
        echo "EDTA PASS LIST IS EMPTY" \\
        > "${modFileName}.EDTA.pass.list"
    fi

    if [ -f "${modFileName}.EDTA.anno/${modFileName}.out" ]; then
        cat "${modFileName}.EDTA.anno/${modFileName}.out" \\
        > "${modFileName}.EDTA.out"
    else
        echo "EDTA DID NOT PRODUCE AN OUT FILE" \\
        > "${modFileName}.EDTA.out"
    fi

    if [ ! -f "${modFileName}.EDTA.TEanno.gff3" ]; then
        echo "##EDTA DID NOT PRODUCE A TEANNO GFF3" \\
        > "${modFileName}.EDTA.TEanno.gff3"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EDTA: \$(EDTA.pl -h | awk ' /##### Extensive/ {print \$7}')
    END_VERSIONS
    """
    
    stub:
    def modFileName = "${fasta_file}.mod"
    """
    touch "${modFileName}.EDTA.TElib.fa"
    touch "${modFileName}.EDTA.intact.gff3"
    touch "${modFileName}.EDTA.pass.list"
    touch "${modFileName}.EDTA.out"
    touch "${modFileName}.EDTA.TEanno.gff3"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EDTA: \$(EDTA.pl -h | awk ' /##### Extensive/ {print \$7}')
    END_VERSIONS
    """
}