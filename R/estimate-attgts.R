

#' Estimate ATT(g,t) point estimates
#' 
#' 
#' @description Given a 'lookup table' of when someone first was treated, calculate 
#' the ATTs of Callaway and Sant'Anna
#'
#' @param g_val Which group to estimate ATT for.
#' @param t_val Which time period to estimate ATT for.
#' @param lookup_table Dataset with "counting" process which reflects when a unit is 
#'      first switched  on.
#' @param N_table Dataset with number of individuals per time and group..
#' @param verbose Whether to return all subcomponents used in ATT calculation (for debugging primarily)
#' @param hetero_var Factor variable for heterogeneous treatment effects. 
#' @param set_div_0_to_0  If there are no individuals present, set probability to 0 instead of NaN 
#'
#' @export 
calculate_att_g_t = function(g_val, 
                             t_val, 
                             lookup_table, 
                             N_table, 
                             hetero_var,
                             verbose = FALSE, 
                             set_div_0_to_0 = FALSE
                             ) {
    if (is.null(hetero_var)) {
        hetero_var = "const"
    }
    if (t_val >= g_val) {
        lag_t_val = g_val - 1
    } else {
        lag_t_val = t_val - 1
    }

    
    type_vals = lookup_table[t == t_val & G == g_val, get(hetero_var)]

    n_g_t_treated = lookup_table[t == t_val & G == g_val, n*w]
    N_g_t = N_table[t == t_val & G == g_val, N*w]
    y_g_t_treated  = n_g_t_treated /  N_g_t

    n_g_gm1_treated = lookup_table[t == lag_t_val & G == g_val, n*w] 
    N_g_gm1 = N_table[t == lag_t_val & G == g_val, N*w]
    y_g_gm1_treated = n_g_gm1_treated / N_g_gm1

    if (all(N_g_gm1 == 0) & set_div_0_to_0 == TRUE) {
        y_g_gm1_treated = 0
    }


    n_g_t_nyt = lookup_table[t == t_val & G != g_val & (t_val < G | G == 0)][, sum(n*w)] 
    N_g_t_nyt = N_table[t == t_val & G != g_val & (t_val < G | G == 0)][, sum(N*w)]
    y_g_t_nyt = n_g_t_nyt / N_g_t_nyt


    n_g_gm1_nyt = lookup_table[t == lag_t_val & G != g_val & (t_val < G | G == 0)][, sum(n*w)] 
    N_g_gm1_nyt = N_table[t == lag_t_val & G != g_val & (t_val < G | G == 0)][, sum(N*w)]
    y_g_gm1_nyt = n_g_gm1_nyt / N_g_gm1_nyt

    if (all(N_g_gm1_nyt == 0) & set_div_0_to_0 == TRUE) {
        y_g_gm1_nyt = 0
    }

    att_g_t = (y_g_t_treated - y_g_gm1_treated) - (y_g_t_nyt - y_g_gm1_nyt)
    if (verbose == TRUE) {
        return(lst(
            n_g_t_treated, 
            N_g_t,
            n_g_gm1_treated,
            N_g_gm1,
            n_g_t_nyt,
            N_g_t_nyt,
            n_g_gm1_nyt, 
            N_g_gm1_nyt,
            att_g_t
            ))
    }
    dt = data.table(
        g = g_val, 
        t = t_val, 
        att_g_t = att_g_t, 
        hetero_var = type_vals
    )
    return(dt)
}