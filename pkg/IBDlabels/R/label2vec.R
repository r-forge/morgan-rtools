label2vec <-
function(label, ngam){

    ##initialise
    cc = label
    flag = 0
    vec = numeric(ngam)
    for( i in (ngam-1):0 ){
        ##f cc = 0 break (flag=0)
        if (cc==0) break 
        ##find gg
        gg = cc - (cc %/% i)*i 
        ##continue if non-zero or cc < CDgg
        if (gg !=0 | cc< maxlabel(i+1) ){ 
            vec[i+1] = gg+1
            cc = cc %/% i
        }else{
            flag =1; break
        }
        i=i-1 
  }
    ##fill rest of s
    if(flag == 0) vec[1:(i+1)] = 1
    if(flag == 1) vec[1:(i+1)] = 1:(i+1)
    return(vec)
}
