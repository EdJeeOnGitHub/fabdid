
rdirichlet = function(n){
    gam_rvar = rgamma(n, 1)
    dir_rvar = gam_rvar / sum(gam_rvar)
    return(dir_rvar)
}

rrademacher = function(n){
    binom_rvar = rbinom(n, 1, 0.5)
    r_rvar = (-1)^binom_rvar
    return(r_rvar)
}

