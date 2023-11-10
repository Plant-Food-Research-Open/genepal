def validateParams(params) {
    validateFastaTags(params)
    validateTETags(params)

    validateTEFastaCorrespondence(params)
}

def validateFastaTags(params) {
    def listOfFastaTuples   = params["target_assemblies"]

    if (isNotListOfLists(listOfFastaTuples, 2)) {
        error 'Error: target_assemblies must be a list of sublists, with each sublist containing 2 elements'
    }

    def fastaTags = listOfFastaTuples.collect { it[0] }

    fastaTags.each {
        if (!(it =~ /^\w+$/)) {
            error "Error: $it is not a valid tag in target_assemblies"
        }
    }

    if (fastaTags.size() != (fastaTags as Set).size()) {
        error "All the tags in target_assemblies should be unique"
    }
}

def validateTETags(params) {
    def listOfTETuples   = params["te_libraries"]

    if (listOfTETuples.isEmpty()) {
        return
    }

    if (isNotListOfLists(listOfTETuples, 2)) {
        error 'Error: te_libraries must be a list of sublists, with each sublist containing 2 elements'
    }

    def teTags = listOfTETuples.collect { it[0] }

    teTags.each {
        if (!(it =~ /^\w+$/)) {
            error "Error: $it is not a valid tag in te_libraries"
        }
    }
}

def validateTEFastaCorrespondence(params) {
    
    def listOfTETuples   = params["te_libraries"]
    def listOfFastaTuples   = params["target_assemblies"]

    def fastaTags = listOfFastaTuples.collect { it[0] }
    def teTags = listOfTETuples.collect { it[0] }

    teTags.each {
        if(!fastaTags.contains(it)) {
            error "Error: $it in te_libraries does not have a corresponding tag in target_assemblies"
        }
    }
}

def isNotListOfLists(thisOne, subListSize) {
    return (!(thisOne instanceof List) || thisOne.isEmpty() || thisOne.any { !(it instanceof List) || it.size() != subListSize })
}