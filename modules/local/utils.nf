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
