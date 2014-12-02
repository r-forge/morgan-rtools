label2jaq <-
function( label, phased = TRUE ){
    
    ##    ngam <- 4
    ##    nlabels <- 15
    
    out <- allJaq()

    jaq <- NULL
    sel <- ifelse( phased, 2,3 )
    for( i in 1:length(label)){

        jaq <- c(jaq, out[ out[,1]==label[i], sel] )
    }
    return( jaq ) 
}
