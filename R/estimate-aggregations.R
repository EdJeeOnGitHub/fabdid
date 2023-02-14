
calculate_event_study = function(att_pr_df){
    manual_es = att_pr_df[, .(wt = pr/sum(pr), att_g_t = att_g_t), event.time][, .(estimate = sum(wt*att_g_t)), event.time ]
    return(manual_es)
}
