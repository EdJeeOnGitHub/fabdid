
create_indiv_first_treat_dt = function(dt, 
                                       y_var, 
                                       group_var, 
                                       id_var,
                                       birth_var = NULL, 
                                       cluster_id = id_var, 
                                       hetero_var = NULL) {
    if (is.null(hetero_var)) {
        hetero_var = "const"
        dt[, const := 1]
    }
    if (!is.null(birth_var)) {
        summ_dt = dt[
            , 
            .(
                first_Y = unique(get(y_var)),
                G = unique(get(group_var)), 
                born_period = unique(get(birth_var)), 
                cluster_id = unique(get(cluster_id)), 
                hetero_var = unique(get(hetero_var))
            ), 
            by = id_var
            ]
        setnames(
            summ_dt, 
            c(id_var, "first_Y", "G", "born_period", "cluster_id", hetero_var)
        )
    } else {
        summ_dt = dt[
            , 
            .(
                first_Y = unique(get(y_var)),
                G = unique(get(group_var)),
                cluster_id = unique(get(cluster_id)),
                hetero_var = unique(get(hetero_var))
            ), 
            by = id_var
            ]
        setnames(
            summ_dt, 
            c(id_var, "first_Y", "G", "cluster_id", hetero_var)
        )

    }
    summ_dt[, rowid := 1:.N]
    return(summ_dt)
}

#' @title Create 'Group' counting data
#'
#' Group individual counting data by treatment group and aggregate by summing.
#'
#' @param dt Indiv "counting" data.
#' @param y_var What is the y variable called.
#' @param group_var What is the group variable called.
#' @param t_levels Unique levels of t to calculate over
#' @param group_levels Unique levels of g to calculate over
#' @param weight_df Weighting df
#' @param hetero_var Factor variable to use to calculate heterogeneous treatment effects 
#'
#' @export
create_group_first_treat_dt = function(dt, 
                                       y_var, 
                                       group_var, 
                                       t_levels, 
                                       group_levels, 
                                       hetero_var = NULL,
                                       weight_df = NULL) {
    
    if (is.null(weight_df)) {
        weight_df = data.frame(
            w = 1,
            G = group_levels
        )
    }
    if (is.null(hetero_var)) {
        hetero_var = "const"
        dt[, const := 1]
    }

    by_vars = c(y_var, group_var, hetero_var)
    summ_group_dt = dt[, .N, by = by_vars][, t := get(y_var)]

    hetero_levels = dt[, unique(get(hetero_var))]

    zero_t_group_dt = CJ(t = t_levels, group_var = group_levels, hetero_var = hetero_levels, n_zero = 0)
    setnames(zero_t_group_dt, c("t", "group_var", hetero_var, "n_zero"))
    full_group_dt = merge(
        summ_group_dt,
        zero_t_group_dt,
        by.x = c("t", group_var, hetero_var),
        by.y = c("t", "group_var", hetero_var), 
        all.y = TRUE
    )
    full_group_dt[is.na(N), N := n_zero]
    full_group_dt[, n_zero := NULL]
    full_cumsum_group_dt = full_group_dt[
        order(get(group_var), t, get(hetero_var)),
        .(n = cumsum(N), t = unique(t)),
        c(group_var, hetero_var)
    ]
    setcolorder(full_cumsum_group_dt, c(group_var, "t", hetero_var, "n"))
    full_cumsum_group_dt = merge(
        full_cumsum_group_dt,
        weight_df,
        by = "G", 
        all.x = TRUE
    )
    return(full_cumsum_group_dt)
}

create_indiv_per_period_dt = function(df, 
                                      group_var, 
                                      t_var, 
                                      t_levels, 
                                      group_levels, 
                                      hetero_var = NULL,
                                      weight_df = NULL) {
    if (is.null(weight_df)) {
        weight_df = data.frame(
            w = 1,
            G = group_levels
        )
    }
    if (is.null(hetero_var)) {
        hetero_var = "const"
        df[, const := 1]
    }

    hetero_levels = df[, unique(get(hetero_var))]

    empty_dt = CJ(G = group_levels, t = t_levels, hetero_var = hetero_levels, N_0 = 0 )
    setnames(empty_dt, c("G", "t", hetero_var, "N_0" ))
    n_indiv = df[, .N, by = c(group_var, t_var, hetero_var)]

    full_dt = merge(
        empty_dt,
        n_indiv,
        all.x = TRUE, 
        by.y = c(group_var, t_var, hetero_var), 
        by.x = c("G", "t", hetero_var)
    )

    full_dt[is.na(N), N := N_0]
    full_dt[, N_0 := NULL]

    full_dt = merge(
        full_dt,
        weight_df,
        by = "G", 
        all.x = TRUE
    )

    pr_df = df[, .(n = .N), by = c(group_var, hetero_var)][, .(pr = n/sum(n), G = get(group_var), hetero_var = get(hetero_var))]
    setnames(
        pr_df, 
        c("pr", "G", hetero_var)
    )
    full_dt = merge(
        full_dt,
        pr_df,
        by = c("G", hetero_var)
    )
    setorderv(full_dt, c("G", "t", hetero_var))
    setkeyv(full_dt, c("G", "t", hetero_var))
    return(full_dt)
}

#' Create N individuals Per Period Using Summary Data
#'
#' @param time_levels time periods to calculate over
#' @param summ_indiv individual level summary data to use
#' @param hetero_var Factor variable for heterogeneous treatment effects. 
#'
#' @export 
create_N_per_period_from_summ = function(time_levels, summ_indiv, hetero_var = NULL){
    if (is.null(hetero_var)) {
        hetero_var = "const"
    }
    N_dt = lapply(
        time_levels,
        function(x){summ_indiv[, .(N = sum(born_period <= x), t = x), c("G", hetero_var)]}
    ) %>%
        rbindlist()
    summ_N_dt = N_dt[, .(N = sum(N)), c("G", hetero_var)][, .(pr = N/sum(N), G, hetero_var = get(hetero_var))]
    setnames(summ_N_dt, c("pr", "G", hetero_var))
    N_dt = merge(
        N_dt, 
        summ_N_dt,
        by.x = c("G", hetero_var), 
        by.y = c("G", hetero_var), 
        all.x = TRUE
    )
    N_dt[, w := 1 ]
    setorderv(N_dt, c("G", "t", hetero_var))
    setcolorder(N_dt, c("G", "t",  "N", "w", "pr"))
    setkeyv(N_dt, c("G", "t", hetero_var))
    return(N_dt)
}