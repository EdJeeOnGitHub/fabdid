
calculate_event_study = function(att_pr_df, y_var){
    manual_es = att_pr_df[!is.na(get(y_var)), .(wt = pr/sum(pr), att_g_t = get(y_var)), event.time][, .(estimate = sum(wt*att_g_t)), event.time ]
    return(manual_es)
}


calculate_group_average = function(att_pr_df, y_var, aggregate_fits = FALSE) {
    if (aggregate_fits == TRUE) {
        group_estim = att_pr_df[
            event.time >= 0 & !is.na(get(y_var))
            ][, 
            .(estimate = sum(get(y_var)*pr_split), pr = sum(pr*pr_split)), 
            group]
    } else {
        group_estim = att_pr_df[event.time >= 0 & !is.na(get(y_var))][, .(estimate = mean(get(y_var)), pr = unique(pr)), group]
    }
    group_estim[, group := as.character(group)]

    ave_group = group_estim[, .(wt = pr/sum(pr), estimate)][, .(group = "Average", estimate = sum(estimate*wt))]
    group_estim[, pr := NULL]
    group_dt = rbindlist(
        list(
            ave_group, 
            group_estim
        )
    )
    return(group_dt)
}
