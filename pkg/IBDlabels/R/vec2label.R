vec2label <-
function(vec, renum = TRUE){

    if( renum ) vec = rplc( vec )
    n = length(vec)
    label = 0
    for( i in 0:(n-1)){
        label = label*i + (vec[i+1] - 1)
    }
    return(label)
}
