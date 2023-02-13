
compare_timings = FALSE

if (compare_timings) {
comp_mb = microbenchmark::microbenchmark(
    # cs_fit = att_gt(
    #     data = binary_sim_df,
    #     yname = "Y_binary",
    #     tname = "period",
    #     gname = "G",
    #     est_method = "ipw",
    #     bstrap = FALSE,
    #     idname = "id",
    #     control_group = "notyettreated" ),
    just_att = map2(
        manual_did$g,
        manual_did$t,
        ~calculate_att_g_t(
            g_val = .x,
            t_val = .y,
            lookup_table = summ_group_dt,
            N_table = N_indiv_dt 
        )
    ),
    manual_did = estimate_did(
        df,
        y_var = "first_Y",
        group_var = "G",
        t_var = "period",
        id_var = "id"
    ), 
    manual_inf = map2(
        manual_did$g, 
        manual_did$t,
        ~calculate_influence_function(
            g_val = .x, 
            t_val = .y, 
            summ_indiv_dt
        )
    ),
    manual_no_probit_inf = map2(
        manual_did$g, 
        manual_did$t,
        ~calculate_influence_function(
            g_val = .x, 
            t_val = .y, 
            summ_indiv_dt, 
            prop_score_known = TRUE
        )
    ),
    times = 2
)
comp_mb
}