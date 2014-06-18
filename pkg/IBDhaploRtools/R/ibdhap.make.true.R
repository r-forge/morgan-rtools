ibdhap.make.true <- function( true.filename){

    ## read in the whole data, remove names and transpose
    rawdat <- read.table( true.filename, sep = " ", header = FALSE )
    outdat <- t( rawdat[,-(1:2)] )
    ## remove the row of NA at the bottom
    dm <- dim( outdat )[1]
    outdat <- outdat[ -dm , ]
    
    return( outdat )
}
