process RESTORE_EDTA_IDS {
    tag "$meta.id"
    label "process_single"

    container "docker://gallvp/python3npkgs:v0.4"

    input:
    tuple val(meta), path(te_lib_fa)
    path(intact_gff3)
    path(pass_list)
    path(out_file)
    path(te_anno_gff3)
    path(renamed_ids_tsv)
    
    output:
    tuple val(meta), path("${meta.id}.EDTA.TElib.fa")               , emit: te_lib_fasta
    tuple val(meta), path("${meta.id}.EDTA.intact.gff3")            , emit: intact_gff3
    tuple val(meta), path("${meta.id}.renamed.ids.EDTA.pass.list")  , emit: pass_list
    tuple val(meta), path("${meta.id}.renamed.ids.EDTA.out")        , emit: out_file
    tuple val(meta), path("${meta.id}.EDTA.TEanno.gff3")            , emit: te_anno_gff3
    tuple val(meta), path("${meta.id}.renamed.ids.tsv")             , emit: renamed_ids_tsv
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $pass_list > "${meta.id}.renamed.ids.EDTA.pass.list"
    cat $out_file > "${meta.id}.renamed.ids.EDTA.out"
    cat $te_lib_fa > "${meta.id}.EDTA.TElib.fa"
    cat $renamed_ids_tsv > "${meta.id}.renamed.ids.tsv"
    
    renamed_ids_head=\$(head -n 1 "$renamed_ids_tsv")
    
    if [[ \$renamed_ids_head == "IDs have acceptable length and character. No change required." ]]; then
        cat $te_anno_gff3 > "${meta.id}.EDTA.TEanno.gff3"
        cat $intact_gff3 > "${meta.id}.EDTA.intact.gff3"
    else
        reverse_edta_naming.py "$renamed_ids_tsv" "$te_anno_gff3" "$intact_gff3" "$meta"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reverse_edta_naming: \$(md5sum \$(which reverse_edta_naming.py) | cut -d' ' -f1)
    END_VERSIONS
    """
    
    stub:
    """
    touch "${meta.id}.EDTA.TElib.fa"
    touch "${meta.id}.EDTA.intact.gff3"
    touch "${meta.id}.renamed.ids.EDTA.pass.list"
    touch "${meta.id}.renamed.ids.EDTA.out"
    touch "${meta.id}.EDTA.TEanno.gff3"
    touch "${meta.id}.renamed.ids.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        reverse_edta_naming: \$(md5sum \$(which reverse_edta_naming.py) | cut -d' ' -f1)
    END_VERSIONS
    """
}