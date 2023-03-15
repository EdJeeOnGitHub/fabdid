

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
#' @param cluster_var Variable to cluster SEs on. 
#' @param weight_df Dataframe of weights.
#' @param prop_score_known Is the propensity score known or estimated
#' @param biter How many bootstrap iterations to use to calculate SEs
#' @param n_cores How many cores to use when bootstrapping
#' @param birth_var If a panel with children being born, what period are they born in.
#' @param hetero_var Factor variable used to estimate heterogeneous treatment effects. 
#' 
#' @export 
estimate_did = function(data, 
                        y_var, 
                        group_var, 
                        t_var, 
                        id_var,
                        birth_var = NULL,
                        cluster_var = NULL,
                        weight_df = NULL,
                        prop_score_known = FALSE, 
                        biter = 1000, 
                        n_cores = 8, 
                        hetero_var = NULL
                        ){
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
        weight_df = weight_df, 
        hetero_var = hetero_var
    )

    if (is.null(cluster_var)) {
        cluster_var = id_var
    }
    summ_indiv_data = create_indiv_first_treat_dt(
        dt = data,
        y_var = y_var,
        group_var = group_var,
        id_var = id_var,
        cluster_id = cluster_var,
        birth_var = birth_var,
        hetero_var = hetero_var
    )
    
    setkeyv(summ_indiv_data, c(group_var, y_var, hetero_var))
    summ_group_data = create_group_first_treat_dt(
        dt = summ_indiv_data,
        y_var = y_var,
        group_var = group_var,
        t_levels = t_levels,
        group_levels = group_levels,
        weight_df = weight_df,
        hetero_var = hetero_var
    )
    setkeyv(summ_group_data, c("G", "t", hetero_var))
    gs_and_ts_we_want = CJ(group = group_levels[group_levels != 0], time = t_levels[t_levels != min(t_levels)])
    # Just put
    att_estimates = map2_dfr(
        gs_and_ts_we_want$group,
        gs_and_ts_we_want$time,
        ~calculate_att_g_t(
            g_val = .x,
            t_val = .y,
            lookup_table = summ_group_data,
            N_table = N_indiv_data, 
            hetero_var = hetero_var
        ), 
        .progress = TRUE
    ) %>%
        as.data.table()

    if (is.null(hetero_var)) {
        hetero_var = "const"
        data[, const := 1]
    }
    gs_ts_heteros_we_want = CJ(
        group = group_levels[group_levels != 0], 
        time = t_levels[t_levels != min(t_levels)], 
        hetero = data[, unique(get(hetero_var))]
        )
    if (is.null(birth_var)) {
        inf_func_output = pmap(
            list(
                gs_ts_heteros_we_want$group,
                gs_ts_heteros_we_want$time,
                gs_ts_heteros_we_want$hetero
            ),
            ~calculate_influence_function(
                g_val = ..1, 
                t_val = ..2, 
                hetero_val = ..3, 
                hetero_var = hetero_var,
                lookup_indiv_table = summ_indiv_data,
                prop_score_known = prop_score_known
            )
        )
    } else {
        inf_func_output = pmap(
            list(
                gs_ts_heteros_we_want$group,
                gs_ts_heteros_we_want$time,
                gs_ts_heteros_we_want$hetero
            ),
            ~calculate_rc_influence_function(
                g_val = ..1,
                t_val = ..2,
                hetero_val = ..3,
                hetero_var = hetero_var,
                lookup_indiv_table = summ_indiv_data,
                row_id_var = "rowid",
                prop_score_known = prop_score_known
            )
        )
    }

    setnames(att_estimates, c("group", "time", "att_g_t", hetero_var))
    att_estimates[, event.time := time - group]
    
    pr_df = N_indiv_data[, .(pr = unique(pr)), c("G", hetero_var)]

    hetero_pr = N_indiv_data[, .(pr = unique(pr)), c("G", hetero_var)][, .(pr_hetero = pr/sum(pr), tmp_het = get(hetero_var)), G]

    att_estimates = merge(
        att_estimates,
        pr_df,
        by.x = c("group", hetero_var),
        by.y = c("G", hetero_var),
        all.x = TRUE
    )
    att_estimates = merge(
        att_estimates,
        hetero_pr,
        by.x = c("group", hetero_var),
        by.y = c("G", "tmp_het"),
        all.x = TRUE
    )
    att_estimates[, treated := event.time >= 0]
    return(lst(att_df = att_estimates, inf_func_output))
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
#' @param cluster_id Nx1 vector of cluster IDs
#' @param hetero_var Factor variable for heterogeneous effects 
#' @export 
estimate_event_study = function(att_df, 
                                inf_matrix,
                                balance_e = NULL, 
                                y_var = "att_g_t", 
                                biter = 1000,
                                n_cores = 8,
                                prop_score_known = FALSE,
                                group_vector = NULL, 
                                cluster_id = NULL, 
                                hetero_var = NULL
                                  ) {
    if (is.null(hetero_var)) {
        indices_we_want = list(1:nrow(att_df))
        hetero_vals = 1
        hetero_var = "const"
        if (prop_score_known == FALSE) {
            stop("Heterogeneous effects and unknown propensity score not yet implemented.")
        }
    } else {
        hetero_vals = att_df[, unique(get(hetero_var))]
        indices_we_want = map(
            hetero_vals, 
            ~which(att_df[, get(hetero_var) == .x])
            )
    }
    es_df = imap_dfr(
        indices_we_want,
        ~.estimate_event_study(
            att_df = att_df[.x, ],
            inf_matrix = inf_matrix[, .x],
            balance_e = balance_e, 
            y_var = y_var,
            biter = biter,
            n_cores = n_cores,
            prop_score_known = prop_score_known,
            group_vector = group_vector, 
            cluster_id = cluster_id
        ) %>%
        as.data.table() %>%
        mutate(tmp_var = hetero_vals[.y])
    )
    setnames(es_df, c("event.time", 
                      "estimate", 
                      "std.error", 
                      "conf.low", 
                      "conf.high", 
                      hetero_var))    
    setorderv(es_df, c("event.time", hetero_var))
    if (hetero_var == "const") {
        es_df[, const := NULL]
    }
    return(es_df)
}
                                


.estimate_event_study = function(att_df, 
                                inf_matrix,
                                balance_e = NULL, 
                                y_var = "att_g_t", 
                                biter = 1000,
                                n_cores = 8,
                                prop_score_known = FALSE,
                                group_vector = NULL, 
                                cluster_id = NULL
                                  ) {
    # att_df = manual_did
    # biter = 100
    # n_cores = 4
    # prop_score_known = FALSE
    # group_vector = summ_indiv_dt[, G]
    # balance_e = NULL
    # y_var = "att_g_t"
    if (!is.null(balance_e)) {
        att_df[, max_et := max(event.time), group]
        att_df = att_df[max_et >= balance_e & event.time <= balance_e]
    }
    # Dealing with NA ATTs
    not_na_atts = att_df[, !is.na(get(y_var))]
    att_df = att_df[not_na_atts]
    inf_matrix = inf_matrix[, not_na_atts]


    att_df[, wt := pr/sum(pr), .(event.time, treated)]
    # agg_pr = att_df[, unique(pr), .(event.time)][, V1]
    agg_pr = att_df[, sum(pr), .(group, event.time)][order(event.time, group), V1]
    es_df = calculate_event_study(att_df, y_var = y_var)

    setorder(es_df, event.time)


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
        weight_if = NULL
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
        ~calculate_se(.x, biter = biter, n_cores = n_cores, cluster_id = cluster_id)
    )
    es_df[, std.error := et_se]
    setorder(es_df, event.time)
    es_df[, conf.low := estimate - 1.96*std.error]
    es_df[, conf.high := estimate + 1.96*std.error]
    return(es_df)
}

#' Estimate Group Average 
#'
#' Given a data.table of ATTs and an influence matrix, create group level estimates.
#' 
#' Assumes propensity score is known.
#' 
#' @param att_df Data.table of ATTs
#' @param inf_matrix Influence function matrix, N individuals by K ATTs
#' @param balance_e Balance event study for e time periods
#' @param y_var Column name in att_df
#' @param biter Bootstrap draws for inference
#' @param n_cores Number of cores for parallal processing
#' @param cluster_id Nx1 vector of cluster IDs
#' @param aggregate_fits  Are multiple fits across covariate factors present in att_df?
#' 
#' @export 
estimate_group_average = function(att_df, 
                                  inf_matrix,
                                  balance_e = NULL, 
                                  y_var = att_g_t, 
                                  biter = 1000,
                                  n_cores = 8,
                                  cluster_id = NULL,
                                  aggregate_fits = FALSE
                                  ) {
    # att_df = manual_did
    # biter = 100
    # n_cores = 4
    # prop_score_known = FALSE
    # group_vector = summ_indiv_dt[, G]
    # balance_e = NULL
    # y_var = "att_g_t"

    # selective.se.inner <- lapply(glist, function(g) {
    #   whichg <- which( (group == g) & (g <= t) & (t<= (group + max_e)))  ### added last condition to allow for limit on longest period included in att
    #   inf.func.g <- as.numeric(get_agg_inf_func(att=att,
    #                                             inffunc1=inffunc1,
    #                                             whichones=whichg,
    #                                             weights.agg=pg[whichg]/sum(pg[whichg]),
    #                                             wif=NULL))
    #   se.g <- getSE(inf.func.g, dp)
    #   list(inf.func=inf.func.g, se=se.g)
    # })
 
    if (!is.null(balance_e)) {
        att_df[, max_et := max(event.time), group]
        att_df = att_df[max_et >= balance_e & event.time <= balance_e]
    }
    group_dt = calculate_group_average(att_df, y_var = y_var, aggregate_fits = aggregate_fits)
    groups = unique(att_df$group)

    if (aggregate_fits == FALSE) {
        pg = att_df[, pr]
        pgg_df = att_df[, .(pr = unique(pr)), .(group, time)]
    } else {
        pg = att_df[, pr*pr_split]
        pgg_df = att_df[, .(pr = sum(pr*pr_split)), .(group, time)]
    }
    # If multiple pr(G = g) per group (i.e. because we've split up the group by 
    # conditioning on a covariate factor, make unique by summing over g&t)
    pgg = pgg_df[match(groups, pgg_df[, group]), pr]

    whichones = map(
        groups,
        ~which(att_df[, group == .x & event.time >= 0])
    )

    inf_func_group = map(
        whichones,
        ~get_agg_inf_func(
            att = att_pr_df,
            inffunc1 = inf_matrix, 
            whichones = .x,
            weights.agg = pg[.x]/sum(pg[.x]) , 
            wif = NULL
        )
    )

    group_se = map_dbl(
        inf_func_group,
        ~calculate_se(
            .x, 
            biter = biter, 
            n_cores = n_cores, 
            cluster_id = cluster_id)
    )

    overall_if_matrix = do.call(cbind, inf_func_group)
    
    # get overall influence function
    overall_if = get_agg_inf_func(att=att_df[group == 'Average', estimate],
                                           inffunc1=overall_if_matrix,
                                           whichones=(1:ncol(overall_if_matrix)),
                                           weights.agg=pgg/sum(pgg),
                                           wif=NULL)


    overall_se = calculate_se(
        overall_if, 

        biter = biter, 
        n_cores = n_cores, 
        cluster_id = cluster_id
    )
    group_dt[, std.error := c(overall_se, group_se)]
    group_dt[, conf.low := estimate - 1.96*std.error]
    group_dt[, conf.high := estimate + 1.96*std.error]
    return(group_dt)
}
