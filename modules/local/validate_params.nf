def validateParams(params) {

    if (!params['repeat_annotator']) {
        error "Error: repeat_annotator must be either 'repeatmodeler' or 'edta'"
    }

    if ( !(params['repeat_annotator'] in ['repeatmodeler', 'edta']) ) {
        error "Error: repeat_annotator must be either 'repeatmodeler' or 'edta'"
    }

    validateRiboDBManifest(params)
    validateLiftoffXrefs(params)
}

def validateRiboDBManifest(params) {
    if (params.remove_ribo_rna) {
        file_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)

        if (file_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${file_ribo_db.getName()}!"}
    }
}

def validateLiftoffXrefs(params) {
    if(!params["liftoff_xref_annotations"]) {
        return
    }

    if(isNotListOfLists(params["liftoff_xref_annotations"]), 2) {
        error "Error: liftoff_xref_annotations must be a list of sublists, with each sublist containing 2 elements"
    }
}

def isNotListOfLists(thisOne, subListSize) {

    if ( subListSize instanceof Integer ) {
        return (!(thisOne instanceof List) || thisOne.isEmpty() || thisOne.any { !(it instanceof List) || it.size() != subListSize })
    }

    return (!(thisOne instanceof List) || thisOne.isEmpty() || thisOne.any { !( it instanceof List ) || !( subListSize.contains( it.size() ) ) })
}

def id_from_file_name(file_name) {

    def trial = ( file_name
        ).replaceFirst(
            /\.f(ast)?q$/, ''
        ).replaceFirst(
            /\.f(asta|sa|a|as|aa)?$/, ''
        ).replaceFirst(
            /\.gff(3)?$/, ''
        ).replaceFirst(
            /\.gz$/, ''
        )

    if ( trial == file_name ) { return file_name }

    return id_from_file_name ( trial )
}
