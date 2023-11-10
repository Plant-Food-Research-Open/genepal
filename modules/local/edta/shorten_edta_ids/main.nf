process SHORTEN_EDTA_IDS {
    tag "$meta.id"
    label "process_single"

    container "docker://gallvp/python3npkgs:v0.4"
    
    input:
    tuple val(meta), path(fasta_file)
    
    output:
    tuple val(meta), path("*.renamed.ids.fa")   , emit: renamed_ids_fasta
    tuple val(meta), path("*.renamed.ids.tsv")  , emit: renamed_ids_tsv
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    FILE="$fasta_file"
    output_prefix="\${FILE%.*}"

    shorten_fasta_ids.py "$fasta_file" "\$output_prefix"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shorten_fasta_ids: \$(md5sum \$(which shorten_fasta_ids.py) | cut -d' ' -f1)
    END_VERSIONS
    """
    
    stub:
    """
    FILE="$fasta_file"
    output_prefix="\${FILE%.*}"

    touch "\${output_prefix}.renamed.ids.fa"
    touch "\${output_prefix}.renamed.ids.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shorten_fasta_ids: \$(md5sum \$(which shorten_fasta_ids.py) | cut -d' ' -f1)
    END_VERSIONS
    """
}