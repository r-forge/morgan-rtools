jaq2label <-
function( jaq, phased = TRUE){

    ##    ngam <- 4 ## only for 4 gametes
    ##    nlabels <- 15

    out <- allJaq()

    label <- NULL
    sel <- ifelse( phased, 2,3 )
    for( i in 1:length(jaq)){

        label <- c( label, out[ out[,sel]==jaq[i], 1] )
    }

    if(!phased){ cat( "WARNING: unphased-jacquard and label are not 1-1 \n" ) }
    
    return( label )
}
