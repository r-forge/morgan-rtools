allLex <-
function( ngam ){
    ## All lexicographic labels

    ## library(mgcv) ## for uniquecombs

    nlabels <- maxlabel( ngam )
    all.vecs <- allVec( ngam )
    all.simp <- t( apply( all.vecs, 1, rplc )  )
    all.lex <- attr( uniquecombs( all.simp ), "index" )

    names( all.lex ) <- 0:nlabels
    return( all.lex )
}
