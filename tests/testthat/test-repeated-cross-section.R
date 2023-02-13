library(BMisc)
library(data.table)
library(dplyr)
library(purrr)
library(furrr)

set.seed(10202)
if (interactive()) {
    devtools::load_all()
}

ncl <- 1
time.periods <- 5
biters <- 200

# Creates simulation params
sim_params = did::reset.sim(
    time.periods = time.periods, 
    n = 2000
)

sim_df = did::build_sim_dataset(
    sp_list = sim_params, 
    panel = TRUE) %>%
    as_tibble()


binary_sim_df = sim_df %>%
    group_by(G) %>%
    mutate(Y_above = Y > quantile(Y, 0.333)) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(
        first_Y =  min(period[Y_above == TRUE]), 
        first_Y = if_else(!is.finite(first_Y), max(period), first_Y)
    ) %>%
    mutate(Y_binary = period >= first_Y)

# introduce births!
binary_sim_df = binary_sim_df %>%
    group_by(id) %>%
    mutate(
        birth_period = round(runif(1, 0, first_Y))
    )  %>%
    ungroup()


rc_sim_df = binary_sim_df %>%
    dplyr::filter(period >= birth_period)



df = as.data.table(rc_sim_df)


tictoc::tic()
cs_fit = did::att_gt(
    data = rc_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    est_method = "ipw",
    idname = "id",
    print_details = FALSE,
    panel = FALSE,
    control_group = "notyettreated" )
tictoc::toc()

tidy_cs_fit = broom::tidy(cs_fit) %>% dplyr::as_tibble()


manual_did = estimate_did(
    df,
    y_var = "first_Y",
    group_var = "G",
    t_var = "period",
    id_var = "id"
)$att_df


comp_df = inner_join(
    manual_did %>% rename(manual_estimate = att_g_t),
    tidy_cs_fit %>%
        select(group, time, cs_estimate = estimate), 
    by = c("group","time")
)



lookup_max_error = comp_df %>%
    dplyr::filter(!(is.na(manual_estimate) & is.na(cs_estimate))) %>%
    mutate(error = abs(manual_estimate - cs_estimate)) %>%
    summarise(max_error = max(error)) %>%
    pull()

test_that("Manual Estimates OK", {
    expect_lte(lookup_max_error, 1e-8)
}
)


# Profiling Stuff
mb_results = microbenchmark::microbenchmark(
    manual_did = estimate_did(data = df, y_var = "first_Y", group_var = 'G', t_var = "period", id_var = "id" ),
    times = 5
)

mean_time_seconds = mean(mb_results$time / 1e9) 

test_that("Timing okay", {
    expect_lte(mean_time_seconds, 5)
})

# Because of nature of birth data we don't actually test against RC here

#### Influence Function Time ####
# N_indiv_dt = create_indiv_per_period_dt(df, "G", "period", c(1:time.periods), unique(df$G))

# summ_indiv_dt = create_indiv_first_treat_dt(
#     dt = df, 
#     y_var = "first_Y", 
#     group_var = "G", 
#     id_var = "id", 
#     birth_var = "birth_period"
# )
# summ_group_dt = create_group_first_treat_dt(
#     summ_indiv_dt, 
#     "first_Y", 
#     "G", 
#     c(1:time.periods), 
#     unique(summ_indiv_dt$G))

# manual_infs = purrr::map2(
#     manual_did$group, 
#     manual_did$time,
#     ~calculate_rc_influence_function(
#         g_val = .x, 
#         t_val = .y, 
#         summ_indiv_dt,
#         row_id_var = "rowid",
#         prop_score_known = TRUE
#     )
# )
# new_infs = map(manual_infs, "full_inf_func")
# new_adjustment = map(manual_infs, "n_adjustment")
# new_infs = map2(new_infs, new_adjustment, ~.x*.y)



# inf_M = matrix(unlist(new_infs), nrow = nrow(summ_indiv_dt))

# wide_comp_ci_df = inner_join(
#     manual_did %>% rename(
#         manual_estimate = att_g_t, 
#         manual_std.error = std.error),
#     tidy_cs_fit %>%
#         select(group, time, cs_estimate = estimate, cs_std.error = std.error), 
#     by = c("group","time")
# )




