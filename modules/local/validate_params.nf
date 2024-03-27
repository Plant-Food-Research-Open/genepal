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

        if (file_ribo_db.isEmpty()) { exit 1, "File provided with --ribo_database_manifest is empty: ${file_ribo_db.getName()}!" }
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

def idFromFileName(fileName) {

    def trial = ( fileName
        ).replaceFirst(
            /\.f(ast)?q$/, ''
        ).replaceFirst(
            /\.f(asta|sa|a|as|aa)?$/, ''
        ).replaceFirst(
            /\.gff(3)?$/, ''
        ).replaceFirst(
            /\.gz$/, ''
        )

    if ( trial == fileName ) { return fileName }

    return idFromFileName ( trial )
}

def validateFastqMetadata(metas, fqs, permAssString) {
    def permAssList = permAssString.split(",")

    // Check if each listed assembly is permissible
    metas.each { meta ->
        if ( meta.target_assemblies.any { !permAssList.contains( it ) } ) {
            exit 1, "Sample ${meta.id} targets ${meta.target_assemblies} which are not in $permAssList"
        }
    }

    // Check if multiple runs of a sample have the same target assemblies
    if ( metas.collect { meta -> meta.target_assemblies }.unique().size() > 1 ) {
        error "Multiple runs of sample ${metas.first().id} must target same assemblies"
    }

    // Check if multiple runs of a sample have the same endedness
    if ( metas.collect { meta -> meta.single_end }.unique().size() > 1 ) {
        error "Multiple runs of sample ${metas.first().id} must have same endedness"
    }

    [ metas.first(), fqs ]
}
