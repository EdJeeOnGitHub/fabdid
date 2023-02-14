library(BMisc)
library(data.table)
library(dplyr)
library(purrr)
library(furrr)

set.seed(289202)

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
        first_Y =  min(period[Y_above == TRUE])
        # first_Y = if_else(!is.finite(first_Y), max(period), first_Y)
    ) %>%
    mutate(Y_binary = period >= first_Y)

# introduce births!
binary_sim_df = binary_sim_df %>%
    group_by(id) %>%
    mutate(
        birth_period = round(runif(1, 0, pmin(first_Y, max(period))))
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
    biters = 1000,
    control_group = "notyettreated" )
tictoc::toc()

tidy_cs_fit = broom::tidy(cs_fit) %>% dplyr::as_tibble()


manual_did = estimate_did(
    df,
    y_var = "first_Y",
    group_var = "G",
    t_var = "period",
    id_var = "id",
    birth_var = "birth_period"
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
    manual_did = estimate_did(
        data = df, 
        y_var = "first_Y", 
        group_var = 'G', 
        t_var = "period", 
        id_var = "id",
        birth_var = "birth_period"
         ),
    times = 5
)

mean_time_seconds = mean(mb_results$time / 1e9) 

test_that("Timing okay", {
    expect_lte(mean_time_seconds, 5)
})

# Because of nature of birth data we don't actually test against RC here

#### Influence Function Time ####
N_indiv_dt = create_indiv_per_period_dt(df, "G", "period", c(1:time.periods), unique(df$G))

summ_indiv_dt = create_indiv_first_treat_dt(
    dt = df, 
    y_var = "first_Y", 
    group_var = "G", 
    id_var = "id", 
    birth_var = "birth_period"
)
summ_group_dt = create_group_first_treat_dt(
    summ_indiv_dt, 
    "first_Y", 
    "G", 
    c(1:time.periods), 
    unique(summ_indiv_dt$G))

manual_infs = purrr::map2(
    manual_did$group, 
    manual_did$time,
    ~calculate_rc_influence_function(
        g_val = .x, 
        t_val = .y, 
        summ_indiv_dt,
        row_id_var = "rowid",
        prop_score_known = FALSE
    )
)


new_infs = map(manual_infs, "full_inf_func")
new_adjustment = map(manual_infs, "n_adjustment")
new_infs = map2(new_infs, new_adjustment, ~.x*.y)

test_that("RC Inf Function Exists", {
    map(
        new_infs,
        ~expect_equal(nrow(.x), nrow(summ_indiv_dt))
    )
})

inf_matrix = matrix(unlist(new_infs), nrow = nrow(summ_indiv_dt))

manual_se = calculate_se(
    inf_matrix,
    biter = 1000
)
panel_cs_fit = did::att_gt(
    data = rc_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    est_method = "ipw",
    idname = "id",
    print_details = FALSE,
    panel = TRUE,
    biters = 1000,
    control_group = "notyettreated" )

cs_panel_se = broom::tidy(panel_cs_fit)$std.error
cs_rc_se = tidy_cs_fit$std.error

comp_se = tibble(
    manual_se = manual_se,
    cs_panel_se = cs_panel_se, 
    cs_rc_se = cs_rc_se
)

comp_se = comp_se %>%
    mutate(
        btw = cs_panel_se < manual_se & manual_se < cs_rc_se
    ) %>%
    mutate(
        bound_diff = if_else(
            btw == TRUE, 
            0, 
            pmin(abs(manual_se - cs_panel_se), abs(manual_se - cs_rc_se))
        )
    ) %>%
    mutate(
        pct_bound_diff = 100*bound_diff/manual_se
    )
test_that("Birth Panel SEs close to True Panel/RC", {
    # We don't stray more than 10% from bound
    expect_lte(mean(comp_se$bound_diff/comp_se$manual_se), 10/100 )
    # And 50% of the time we're within RC and Panel SEs
    expect_gte(mean(comp_se$btw), 0.5)
})

