
create_indiv_first_treat_dt = function(dt, 
                                       y_var, 
                                       group_var, 
                                       id_var,
                                       birth_var = NULL) {
    if (!is.null(birth_var)) {
        summ_dt = dt[
            , 
            .(
                first_Y = unique(get(y_var)),
                G = unique(get(group_var)), 
                born_period = unique(get(birth_var))
            ), 
            by = id_var
            ]
    } else {
        summ_dt = dt[
            , 
            .(
                first_Y = unique(get(y_var)),
                G = unique(get(group_var))
            ), 
            by = id_var
            ]

    }
    summ_dt[, rowid := 1:.N]
    return(summ_dt)
}


create_group_first_treat_dt = function(dt, y_var, group_var, t_levels, group_levels, weight_df = NULL) {
    
    if (is.null(weight_df)) {
        weight_df = data.frame(
            w = 1,
            G = group_levels
        )
    }


    summ_group_dt = dt[, .N, by = c(y_var, group_var)][, t := get(y_var)]
     
    
    zero_t_group_dt = CJ(t = t_levels, group_var = group_levels, n_zero = 0)
    full_group_dt = merge(
        summ_group_dt,
        zero_t_group_dt,
        by.x = c("t", group_var),
        by.y = c("t", "group_var"), 
        all.y = TRUE
    )
    full_group_dt[is.na(N), N := n_zero]
    full_group_dt[, n_zero := NULL]
    full_cumsum_group_dt = full_group_dt[
        order(get(group_var), t),
        .(n = cumsum(N), t = unique(t)),
        group_var
    ]
    setcolorder(full_cumsum_group_dt, c(group_var, "t", "n"))
    full_cumsum_group_dt[, n_lag := shift(n), by = G]
    full_cumsum_group_dt = merge(
        full_cumsum_group_dt,
        weight_df,
        by = "G", 
        all.x = TRUE
    )
    return(full_cumsum_group_dt)
}

create_indiv_per_period_dt = function(df, group_var, t_var, t_levels, group_levels, weight_df = NULL) {
    if (is.null(weight_df)) {
        weight_df = data.frame(
            w = 1,
            G = group_levels
        )
    }
    empty_dt = CJ(G = group_levels, t = t_levels, N_0 = 0 )
    n_indiv = df[, .N, by = c(group_var, t_var)]

    full_dt = merge(
        empty_dt,
        n_indiv,
        all.x = TRUE, 
        by.y = c(group_var, t_var), 
        by.x = c("G", "t")
    )

    full_dt[is.na(N), N := N_0]
    full_dt[, N_0 := NULL]
    full_dt[, N_lag := shift(N), by = group_var]

    full_dt = merge(
        full_dt,
        weight_df,
        by = "G", 
        all.x = TRUE
    )
    pr_df = df[, .(n = .N), by = group_var][, .(pr = n/sum(n), G = get(group_var))]
    full_dt = merge(
        full_dt,
        pr_df,
        by = "G"
    )
    return(full_dt)
}