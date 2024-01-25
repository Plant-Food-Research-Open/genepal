process EVIGENE_TR2AACDS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/evigene:23.7.15--hdfd78af_1':
        'biocontainers/evigene:23.7.15--hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("dropset/*")  , emit: dropset
    tuple val(meta), path("okayset/*")  , emit: okayset
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args     ?: ''
    def prefix      = task.ext.prefix   ?: "${meta.id}"

    def max_memory  = 7*1024
    if (!task.memory) {
        log.info '[evigene] Available memory not known - defaulting to 7GB. Specify process memory requirements to change this.'
    } else {
        max_memory  = (task.memory.mega*0.8).intValue()
    }

    def simple_name = fasta.simpleName
    """
    /usr/local/opt/evigene-23.7.15/scripts/prot/tr2aacds.pl \\
        $args \\
        -NCPU=$task.cpus \\
        -MAXMEM=$max_memory \\
        -cdnaseq $fasta

    find \\
        dropset \\
        -type f \\
        -exec sh -c 'mv "\$1" "\$(echo \$1 | sed s/$simple_name/$prefix/)"' sh {} \\;

    find \\
        okayset \\
        -type f \\
        -exec sh -c 'mv "\$1" "\$(echo \$1 | sed s/$simple_name/$prefix/)"' sh {} \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tr2aacds: \$(cat /usr/local/opt/evigene-23.7.15/scripts/prot/tr2aacds.pl | grep 'use constant VERSION =>  ' | sed 's/use constant VERSION =>  //' | sed 's/;.*//' | sed "s/'//g")
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def max_memory  = 7*1024
    if (!task.memory) {
        log.info '[evigene] Available memory not known - defaulting to 7GB. Specify process memory requirements to change this.'
    }
    """
    mkdir dropset
    touch dropset/${prefix}.drop.aa
    touch dropset/${prefix}.drop.cds
    touch dropset/${prefix}.drop.tr

    mkdir okayset
    touch okayset/${prefix}.ann.txt
    touch okayset/${prefix}.cull.aa
    touch okayset/${prefix}.cull.cds
    touch okayset/${prefix}.cull.mrna
    touch okayset/${prefix}.genesum.txt
    touch okayset/${prefix}.mainalt.tab
    touch okayset/${prefix}.okay.aa
    touch okayset/${prefix}.okay.cds
    touch okayset/${prefix}.okay.mrna
    touch okayset/${prefix}.pubids
    touch okayset/${prefix}.pubids.old
    touch okayset/${prefix}.pubids.realt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tr2aacds: \$(cat /usr/local/opt/evigene-23.7.15/scripts/prot/tr2aacds.pl | grep 'use constant VERSION =>  ' | sed 's/use constant VERSION =>  //' | sed 's/;.*//' | sed "s/'//g")
    END_VERSIONS
    """
}
