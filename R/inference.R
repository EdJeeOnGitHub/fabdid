
#' @title Compute extra term in influence function due to estimating weights
#'
#' @description A function to compute the extra term that shows up in the
#'  influence function for aggregated treatment effect parameters
#'  due to estimating the weights.
#'  TAKEN DIRECTLY FROM DID PACKAGE
#' 
#' @param keepers a vector of indices for which group-time average
#'  treatment effects are used to compute a particular aggregated parameter
#' @param pg a vector with same length as total number of group-time average
#'  treatment effects that contains the probability of being in particular group
#' @param weights.ind additional sampling weights (nx1)
#' @param G vector containing which group a unit belongs to (nx1)
#' @param group vector of groups
#'
#' @return nxk influence function matrix
#'
wif <- function(keepers, pg, weights.ind, G, group) {
  # note: weights are all of the form P(G=g|cond)/sum_cond(P(G=g|cond))
  # this is equal to P(G=g)/sum_cond(P(G=g)) which simplifies things here

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (weights.ind * 1*BMisc::TorF(G==group[k]) - pg[k]) /
      sum(pg[keepers])
  })
  # effect of estimating weights in the denominator
  if2 <- base::rowSums( sapply( keepers, function(k) {
    weights.ind*1*BMisc::TorF(G==group[k]) - pg[k]
  })) %*%
    t(pg[keepers]/(sum(pg[keepers])^2))

  # return the influence function for the weights
  if1 - if2
}

#' Calculate 'rc' style influence function
#'
#'
#' @param g_val Group in ATT(g, t)
#' @param t_val t in ATT(g, t)
#' @param lookup_indiv_table "counting" binary data
#' @param row_id_var Variable tracking individual ID, must be contiguous hence 'rowid'
#' @param verbose Whether to return intermediate outputs
#' @param check Internal debugging function
#' @param prop_score_known Is the propensity score known
#' @param hetero_var Factor variable for heterogeneous treatment effects. 
#' @param hetero_val Level of factor variable to estimate IF for. 
#' 
#' @export
calculate_rc_influence_function = function(g_val, 
                                           t_val, 
                                           hetero_val = NULL,
                                           hetero_var = NULL,
                                           lookup_indiv_table,
                                           row_id_var,
                                           verbose = FALSE,
                                           check = FALSE,
                                           prop_score_known = FALSE) {

    # g_val = 5
    # t_val = 5
    # lookup_table = summ_group_dt
    # N_table = N_indiv_dt
    # lookup_indiv_table = summ_indiv_dt
    # row_id_var = "rowid"
    # prop_score_known = TRUE
    if (t_val >= g_val) {
        lag_t_val = g_val - 1
    } else {
        lag_t_val = t_val - 1
    }

    if (is.null(hetero_var) & is.null(hetero_val)) {
      hetero_var = "const"
      lookup_indiv_table[, const := 1]
      hetero_val = 1
    }


    people_we_want = lookup_indiv_table[, ((G == g_val & get(hetero_var) == hetero_val) | (t_val < G | G == 0)) & born_period <= t_val]
    subset_lookup_indiv_table = lookup_indiv_table[people_we_want]
    subset_lookup_indiv_table[, treated := factor(G == g_val, levels = c(TRUE, FALSE))]
    subset_lookup_indiv_table[, Y_post := first_Y <= t_val]
    subset_lookup_indiv_table[, Y_pre := first_Y <= lag_t_val]
    deltaY = subset_lookup_indiv_table[, Y_post - Y_pre]


    Y_post = subset_lookup_indiv_table[, Y_post]
    Y_pre = subset_lookup_indiv_table[born_period < g_val, Y_pre]

    rc_ids = c(
      subset_lookup_indiv_table[, get(row_id_var)], 
      subset_lookup_indiv_table[born_period < g_val, get(row_id_var)])
    y = c(
        Y_post,
        Y_pre
    )

    post = c(rep(1, length(Y_post)), rep(0, length(Y_pre)))

    D = c(
        subset_lookup_indiv_table[, as.logical(treated)],
        subset_lookup_indiv_table[born_period < g_val, as.logical(treated)]
        )
    n = length(D)
    # if empty group return NA
    if (n == 0) {
      return(lst(g = g_val, t = t_val, NA, n_adjustment = NA))
    }

    i.weights = NULL
    int.cov = matrix(1, n)


  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")

  if (!prop_score_known) {
    #Pscore estimation (logit) and also its fitted values
    PS <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights))
    ps.fit <- as.vector(PS$fitted.values)
    # Do not divide by zero
    ps.fit <- pmin(ps.fit, 1 - 1e-16)
    w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
    w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)
  }
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  #Compute IPW estimator
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  if (prop_score_known) {
    w.cont.pre <- i.weights * (1 - D) * (1 - post)
    w.cont.post <- i.weights * (1 - D) * post
  }

  mu_w.treat.pre = mean(w.treat.pre) + 1e-8
  mu_w.treat.post = mean(w.treat.post) + 1e-8

  mu_w.cont.pre = mean(w.cont.pre) + 1e-8
  mu_w.cont.post = mean(w.cont.post) + 1e-8


  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * y / mu_w.treat.pre
  eta.treat.post <- w.treat.post * y / mu_w.treat.post
  eta.cont.pre <- w.cont.pre * y / mu_w.cont.pre
  eta.cont.post <- w.cont.post * y / mu_w.cont.post

  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)

  # ATT estimator
  ipw.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  if (!prop_score_known) {
    # Asymptotic linear representation of logit's beta's
    score.ps <- i.weights * (D - ps.fit) * int.cov
    Hessian.ps <- stats::vcov(PS) * n
    asy.lin.rep.ps <-  score.ps %*% Hessian.ps
    # Estimation effect from gamma hat (pscore)
    # Derivative matrix (k x 1 vector)
    M2.pre <- base::colMeans(w.cont.pre *(y - att.cont.pre) * int.cov)/mean(w.cont.pre)
    M2.post <- base::colMeans(w.cont.post *(y - att.cont.post) * int.cov)/mean(w.cont.post)
    # Now the influence function related to estimation effect of pscores
    inf.cont.ps <- asy.lin.rep.ps %*% (M2.post - M2.pre)

  } else {
    inf.cont.ps = 0
  }
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.pre <- eta.treat.pre - w.treat.pre * att.treat.pre/ mu_w.treat.pre
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post/mu_w.treat.post
  inf.treat <- inf.treat.post - inf.treat.pre
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf.cont.pre <- eta.cont.pre - w.cont.pre * att.cont.pre/mu_w.cont.pre
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post/mu_w.cont.post
  inf.cont <- inf.cont.post - inf.cont.pre


  # Influence function for the control component
  inf.cont <- inf.cont + inf.cont.ps

  #get the influence function of the DR estimator (put all pieces together)
  att.inf.func <- inf.treat - inf.cont



    full_inf_func = matrix(0, nrow(lookup_indiv_table))
    # We've calculated IF for each "section" separately
    # now mash them together to get each indiv's contribution
    indiv_level_inf_func = aggregate(att.inf.func, list(rc_ids), sum)
    full_inf_func[indiv_level_inf_func[, 1]] = indiv_level_inf_func[, 2] 


    n_all = nrow(lookup_indiv_table)
    n_subset = nrow(subset_lookup_indiv_table)
    return(lst(g = g_val, t = t_val, full_inf_func, n_adjustment = n_all/length(att.inf.func)))
}
test_calculate_influence_function = function(g_val, 
                                        t_val, 
                                        hetero_val = NULL,
                                        hetero_var = NULL,
                                        lookup_indiv_table,
                                        verbose = FALSE,
                                        check = FALSE,
                                        prop_score_known = FALSE) {

    browser()
    # g_val = 3
    # t_val = 3
    # lookup_table = summ_group_dt
    # N_table = N_indiv_dt
    # lookup_indiv_table = summ_indiv_dt
    if (t_val >= g_val) {
        lag_t_val = g_val - 1
    } else {
        lag_t_val = t_val - 1
    }
    
    if (is.null(hetero_var) & is.null(hetero_val)) {
      hetero_var = "const"
      lookup_indiv_table[, const := 1]
      hetero_val = 1
    }


    people_we_want = lookup_indiv_table[, (G == g_val & get(hetero_var) == hetero_val) | (t_val < G | G == 0)]
    subset_lookup_indiv_table = lookup_indiv_table[(G == g_val & get(hetero_var) == hetero_val) | (t_val < G | G == 0)]
    subset_lookup_indiv_table[, treated := factor(G == g_val, levels = c(TRUE, FALSE))]
    pr_treat = subset_lookup_indiv_table[, mean(treated == TRUE)]
    subset_lookup_indiv_table[, Y_post := first_Y <= t_val]
    subset_lookup_indiv_table[, Y_pre := first_Y <= lag_t_val]
    deltaY = subset_lookup_indiv_table[, Y_post - Y_pre]




    n_all =  nrow(lookup_indiv_table)
    n_subset = nrow(subset_lookup_indiv_table)

    D = subset_lookup_indiv_table[, as.logical(treated)]



    n = nrow(subset_lookup_indiv_table)
    if (prop_score_known == FALSE){
        PS = stats::glm(D ~ 1, family = "binomial")
        ps.fit = as.vector(PS$fitted.value)
        w.cont = ps.fit * (1 - D) / (1 - ps.fit)    
    } else {
        w.cont = (1 - D)
    }

    w.treat = D



    att.treat = w.treat*deltaY
    att.cont = w.cont*deltaY


    eta.treat = mean(att.treat) / mean(w.treat)
    eta.cont = mean(att.cont) / mean(w.cont)

    inf.treat = (att.treat - w.treat * eta.treat) / mean(w.treat)
    inf.cont.1 = (att.cont - w.cont*eta.cont)

    if (prop_score_known == FALSE) {
        score.ps = (D - ps.fit)
        Hessian.ps = stats::vcov(PS) * n
        asy.lin.rep.ps = score.ps %*% Hessian.ps
        M2 = mean(w.cont * (deltaY - eta.cont))
        inf.cont.2 = asy.lin.rep.ps %*% M2
    } else {
        inf.cont.2 = matrix(0, n_subset)
    }

    inf.control = (inf.cont.1 + inf.cont.2) / mean(w.cont)
    att.inf.func = inf.treat - inf.control
    att.inf.func = att.inf.func[, 1]

    rowids = which(people_we_want == TRUE, arr.ind = FALSE, useNames = TRUE)
    names(att.inf.func) = rowids

    full_inf_func = matrix(0, nrow(lookup_indiv_table))
    full_inf_func[rowids] = att.inf.func
    return(lst(g = g_val, t = t_val, full_inf_func, n_adjustment = n_all/n_subset))
}


#' Calculate 'rc' style influence function
#'
#' @inherit calculate_rc_influence_function
#' 
#' @export
calculate_influence_function = function(g_val, 
                                        t_val, 
                                        hetero_val = NULL,
                                        hetero_var = NULL,
                                        lookup_indiv_table,
                                        verbose = FALSE,
                                        check = FALSE,
                                        prop_score_known = FALSE) {

    # g_val = 3
    # t_val = 3
    # lookup_table = summ_group_dt
    # N_table = N_indiv_dt
    # lookup_indiv_table = summ_indiv_dt
    if (t_val >= g_val) {
        lag_t_val = g_val - 1
    } else {
        lag_t_val = t_val - 1
    }
    
    if (is.null(hetero_var) & is.null(hetero_val)) {
      hetero_var = "const"
      lookup_indiv_table[, const := 1]
      hetero_val = 1
    }


    people_we_want = lookup_indiv_table[, (G == g_val & get(hetero_var) == hetero_val) | (t_val < G | G == 0)]
    subset_lookup_indiv_table = lookup_indiv_table[(G == g_val & get(hetero_var) == hetero_val) | (t_val < G | G == 0)]
    subset_lookup_indiv_table[, treated := factor(G == g_val, levels = c(TRUE, FALSE))]
    pr_treat = subset_lookup_indiv_table[, mean(treated == TRUE)]
    subset_lookup_indiv_table[, Y_post := first_Y <= t_val]
    subset_lookup_indiv_table[, Y_pre := first_Y <= lag_t_val]
    deltaY = subset_lookup_indiv_table[, Y_post - Y_pre]




    n_all =  nrow(lookup_indiv_table)
    n_subset = nrow(subset_lookup_indiv_table)

    D = subset_lookup_indiv_table[, as.logical(treated)]



    n = nrow(subset_lookup_indiv_table)
    if (prop_score_known == FALSE){
        PS = stats::glm(D ~ 1, family = "binomial")
        ps.fit = as.vector(PS$fitted.value)
        w.cont = ps.fit * (1 - D) / (1 - ps.fit)    
    } else {
        w.cont = (1 - D)
    }

    w.treat = D



    att.treat = w.treat*deltaY
    att.cont = w.cont*deltaY


    eta.treat = mean(att.treat) / mean(w.treat)
    eta.cont = mean(att.cont) / mean(w.cont)

    inf.treat = (att.treat - w.treat * eta.treat) / mean(w.treat)
    inf.cont.1 = (att.cont - w.cont*eta.cont)

    if (prop_score_known == FALSE) {
        score.ps = (D - ps.fit)
        Hessian.ps = stats::vcov(PS) * n
        asy.lin.rep.ps = score.ps %*% Hessian.ps
        M2 = mean(w.cont * (deltaY - eta.cont))
        inf.cont.2 = asy.lin.rep.ps %*% M2
    } else {
        inf.cont.2 = matrix(0, n_subset)
    }

    inf.control = (inf.cont.1 + inf.cont.2) / mean(w.cont)
    att.inf.func = inf.treat - inf.control
    att.inf.func = att.inf.func[, 1]

    rowids = which(people_we_want == TRUE, arr.ind = FALSE, useNames = TRUE)
    names(att.inf.func) = rowids

    full_inf_func = matrix(0, nrow(lookup_indiv_table))
    full_inf_func[rowids] = att.inf.func
    return(lst(g = g_val, t = t_val, full_inf_func, n_adjustment = n_all/n_subset))
}

#' Calculate Standard Errors Using Multiplier Bootstrap
#' 
#' @description Calculate standard errors using influence function
#' 
#' @param inf_matrix NxK matrix with N indiv observations and K ATTs
#' @param biter Number of bootstrap iterations
#' @param pl Run in parallel
#' @param n_cores Number of cores to use in parallel
#' @param alp Test size, defaults to 0.05.
#' @param cluster_id Vector nx1 of cluster IDs. 
#' @export
calculate_se = function(inf_matrix, cluster_id = NULL, biter = 2000, pl = TRUE, n_cores = 8, alp = 0.05){
    if (is.null(cluster_id)) {
        n = nrow(inf_matrix)
        bres = sqrt(n) * run_multiplier_bootstrap(
            inf_matrix, 
            biter, 
            pl = pl, 
            n_cores)
        n_clusters = n
    } else {
        n_clusters = length(unique(cluster_id))
        cluster_n = aggregate(cluster_id, by = list(cluster_id), length)[, 2]
        cluster_mean_if = rowsum(inf_matrix, cluster_id, reorder = TRUE) / cluster_n
        bres = sqrt(n_clusters) * run_multiplier_bootstrap(cluster_mean_if, biter, pl = pl, n_cores)
    }

        V = cov(bres)
        bSigma <- apply(bres, 2,
                        function(b) (quantile(b, .75, type=1, na.rm = T) -
                                        quantile(b, .25, type=1, na.rm = T))/(qnorm(.75) - qnorm(.25)))
        # critical value for uniform confidence band
        bT <- base::suppressWarnings(apply(bres, 1, function(b) max( abs(b/bSigma), na.rm = T)))
        bT <- bT[is.finite(bT)]
        crit.val <- quantile(bT, 1-alp, type=1, na.rm = T)
        se = as.numeric(bSigma) / sqrt(n_clusters)
    return(se)
}

#' Run Multiplier Bootstrap
#'
#' @description Multiplier Bootstrap
#'
#' @param inf.func Influence function matrix. NxK where N is N indiv and K is # of ATTs
#' @param biters Bootstrap iterations. N.B. plural use of arg here but biter elsewhere. Nice. 
#' @param pl Process in parallel
#' @param cores Number of cores to use if parallel processing
#' 
#' TAKEN DIRECTLY FROM BCALLAWAY11/DID
#' @export
run_multiplier_bootstrap <- function(inf.func, biters, pl = FALSE, cores = 1) {
  ngroups = ceiling(biters/cores)
  chunks = rep(ngroups, cores)
  # Round down so you end up with the right number of biters
  chunks[1] = chunks[1] + biters - sum(chunks)

  n <- nrow(inf.func)
  parallel.function <- function(biters) {
    BMisc::multiplier_bootstrap(inf.func, biters)
  }
  # From tests, this is about where it becomes worth it to parallelize
  if(n > 2500 & pl == TRUE & cores > 1) {
    results = parallel::mclapply(
      chunks,
      FUN = parallel.function,
      mc.cores = cores
    )
    results = do.call(rbind, results)
  } else {
    results = BMisc::multiplier_bootstrap(inf.func, biters)
  }
  return(results)
}

#' @title Get an influence function for particular aggregate parameters
#'
#' @description This is a generic internal function for combining influence
#'  functions across ATT(g,t)'s to return an influence function for
#'  various aggregated treatment effect parameters.
#' 
#'
#' @param att vector of group-time average treatment effects
#' @param inffunc1 influence function for all group-time average treatment effects
#'  (matrix)
#' @param whichones which elements of att will be used to compute the aggregated
#'  treatment effect parameter
#' @param weights.agg the weights to apply to each element of att(whichones);
#'  should have the same dimension as att(whichones)
#' @param wif extra influence function term coming from estimating the weights;
#'  should be n x k matrix where k is dimension of whichones
#'
#'  TAKEN DIRECTLY FROM BCALLAWAY11/DID
#' 
#' @return nx1 influence function
#' @export
get_agg_inf_func <- function(att, inffunc1, whichones, weights.agg, wif=NULL) {
  # enforce weights are in matrix form
  weights.agg <- as.matrix(weights.agg)

  # multiplies influence function times weights and sums to get vector of weighted IF (of length n)
  thisinffunc <- inffunc1[,whichones]%*%weights.agg

  # Incorporate influence function of the weights
  if (!is.null(wif)) {
    thisinffunc <- thisinffunc + wif%*%as.matrix(att[whichones])
  }

  # return influence function
  return(thisinffunc)
}
