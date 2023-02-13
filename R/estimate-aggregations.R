
calculate_event_study = function(att_pr_df){
    manual_es = att_pr_df %>%
        group_by(
            event.time
        ) %>%
        mutate(
            wt = pr/sum(pr)
        ) %>%
        summarise(
            estimate = sum(wt*att_g_t)
        )
    return(manual_es)
}
