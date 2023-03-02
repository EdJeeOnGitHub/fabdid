
calculate_event_study = function(att_pr_df, y_var){
    manual_es = att_pr_df[, .(wt = pr/sum(pr), att_g_t = get(y_var)), event.time][, .(estimate = sum(wt*att_g_t, na.rm = TRUE)), event.time ]
    return(manual_es)
}
