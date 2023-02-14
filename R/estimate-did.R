

#' Estimate CS DiD with Binary, Absorbing Treatment
#'
#' Estimates CS DiD when Y data is binary and 'absorbing' (once you're vaccinated
#' you cannot become unvaccinated).
#'
#' @param data Dataset to use for estimation
#' @param y_var Integer Y variable that denotes the *first* time period an individual is vaccinated
#' @param group_var Treatment group, in the style of CS package, must be first time period group is treated
#' @param t_var Time period variable.
#' @param id_var Unique ID for individuals.
#' @param weight_df Dataframe of weights.
#' @param prop_score_known Is the propensity score known or estimated
#' @param biter How many bootstrap iterations to use to calculate SEs
#' @param n_cores How many cores to use when bootstrapping
#' @param birth_var If a panel with children being born, what period are they born in.
#' 
#' @export 
estimate_did = function(data, 
                        y_var, 
                        group_var, 
                        t_var, 
                        id_var,
                        birth_var = NULL,
                        weight_df = NULL,
                        prop_score_known = FALSE, 
                        biter = 1000, 
                        n_cores = 8){
    data = as.data.table(data)
    t_levels = sort(unique(data[, get(t_var)]))
    group_levels = data[,.(G = unique(get(group_var)))][order(G), G]
    setkeyv(data, c(id_var, group_var))
    N_indiv_data = create_indiv_per_period_dt(
        df = data,
        group_var = group_var,
        t_var = t_var,
        t_levels = t_levels, 
        group_levels = group_levels,
        weight_df = weight_df
    )

    summ_indiv_data = create_indiv_first_treat_dt(
        dt = data,
        y_var = y_var,
        group_var = group_var,
        id_var = id_var,
        birth_var = birth_var
    )
    
    setkeyv(summ_indiv_data, c(group_var, y_var))
    summ_group_data = create_group_first_treat_dt(
        dt = summ_indiv_data,
        y_var = y_var,
        group_var = group_var,
        t_levels = t_levels,
        group_levels = group_levels,
        weight_df = weight_df
    )
    setkeyv(summ_group_data, c("G", "t"))
    gs_and_ts_we_want = CJ(group = group_levels[group_levels != 0], time = t_levels[t_levels != min(t_levels)])
    # Just put
    att_estimates = map2(
        gs_and_ts_we_want$group,
        gs_and_ts_we_want$time,
        ~calculate_att_g_t(
            g_val = .x,
            t_val = .y,
            lookup_table = summ_group_data,
            N_table = N_indiv_data
        ), 
        .progress = TRUE
    )

    if (is.null(birth_var)) {
        inf_func_output = map2(
            gs_and_ts_we_want$g, 
            gs_and_ts_we_want$t,
            ~calculate_influence_function(
                g_val = .x, 
                t_val = .y, 
                summ_indiv_data,
                prop_score_known = prop_score_known
            )
        )
    } else {
        inf_func_output = map2(
            gs_and_ts_we_want$g, 
            gs_and_ts_we_want$t, 
            ~calculate_rc_influence_function(
                g_val = .x,
                t_val = .y,
                summ_indiv_data,
                row_id_var = "rowid",
                prop_score_known = prop_score_known
            )
        )
    }


    gs_and_ts_we_want[, att_g_t := map_dbl(att_estimates, "att_g_t")]
    gs_and_ts_we_want[, event.time := time - group]
    
    pr_df = N_indiv_data[, .(pr = unique(pr)), G]

    gs_and_ts_we_want = merge(
        gs_and_ts_we_want,
        pr_df,
        by.x = "group",
        by.y = "G",
        all.x = TRUE
    )
    gs_and_ts_we_want[, treated := event.time >= 0]
    return(lst(att_df = gs_and_ts_we_want, inf_func_output))
}



#' Estimate Event Study
#'
#' Given a data.table of ATTs and an influence matrix, create event study estimates.
#' 
#' @param att_df Data.table of ATTs
#' @param inf_matrix Influence function matrix, N individuals by K ATTs
#' @param balance_e Balance event study for e time periods
#' @param y_var Column name in att_df
#' @param biter Bootstrap draws for inference
#' @param n_cores Number of cores for parallal processing
#' @param prop_score_known Propensity score known or to be estimated
#' @param group_vector Nx1 vector of group IDs

#' @export 
estimate_event_study = function(att_df, 
                                  inf_matrix,
                                  balance_e = NULL, 
                                  y_var = att_g_t, 
                                  biter = 1000,
                                  n_cores = 8,
                                  prop_score_known = FALSE,
                                  group_vector = NULL
                                  ) {
    # att_df = manual_did
    # biter = 100
    # n_cores = 4
    # prop_score_known = FALSE
    # group_vector = summ_indiv_dt[, G]
    if (!is.null(balance_e)) {
        att_df = att_df %>%
            group_by(group) %>%
            mutate(max_et = mean(event.time)) %>%
            filter(max_et >= balance_e) %>%
            filter(event.time <= balance_e) %>%
            ungroup() %>%
            as.data.table()
    }
    att_df[, wt := pr/sum(pr), .(event.time, treated)]
    agg_pr = att_df[, unique(pr), .(event.time)][, V1]
    es_df = att_df[, .(estimate = sum(get(y_var)*wt)), .(event.time, treated)]

    event_times = unique(es_df$event.time)
    event_time_att_idx = map(
        event_times, 
        ~lst(
            whichones = which(att_df[, event.time == .x]), 
            weights.agg = att_df[event.time == .x, wt],
            pg = att_df[event.time == .x, pr]
            )
        )
    if (prop_score_known == FALSE) {
        weight_if = map(
            event_time_att_idx,
            ~wif(
                keepers = .x$whichones,
                pg = agg_pr,
                weights.ind = rep(1, nrow(inf_matrix)),
                G = group_vector,
                group = att_df[, group]
            )
        )
    } else {
        weight_if = map(event_time_att_idx, NULL)
    }


    et_if = imap(
        event_time_att_idx, 
        ~get_agg_inf_func(
            att = att_df[, get(y_var)], 
            inffunc1 = inf_matrix,
            whichones = .x$whichones,
            weights.agg = .x$weights.agg,
            wif = weight_if[[.y]]
        )
    )
    et_se = map_dbl(
        et_if,
        ~calculate_se(.x, biter = biter, n_cores = n_cores)
    )
    es_df[, std.error := et_se]
    setorder(es_df, event.time)
    return(es_df)
}