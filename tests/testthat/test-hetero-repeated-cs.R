
library(BMisc)
library(data.table)
library(dplyr)
library(purrr)
library(furrr)
library(did)

set.seed(289202)

custom_min = function(x) {if (length(x) > 0) min(x) else Inf}
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


probs = c(1:10)/sum(1:10)

binary_sim_df = sim_df %>%
    group_by(G) %>%
    mutate(Y_above = Y > quantile(Y, 0.333)) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(
        first_Y =  custom_min(period[Y_above == TRUE])
        # first_Y = if_else(!is.finite(first_Y), max(period), first_Y)
    ) %>%
    mutate(Y_binary = period >= first_Y) %>%
    group_by(id) %>%
    mutate(
        type = sample(1:10, size = 1, prob = probs)
    )

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


het_manual_fit = estimate_did(
    df,
    y_var = "first_Y",
    group_var = "G",
    t_var = "period",
    id_var = "id",
    birth_var = "birth_period",
    hetero_var = "type"
)$att_df

# Because of nature of birth data we don't actually test against RC here

#### Influence Function Time ####
N_indiv_dt = create_indiv_per_period_dt(
    df, 
    "G", 
    "period", 
    c(1:time.periods), 
    unique(df$G),
    hetero_var = "type"
    )


summ_indiv_dt = create_indiv_first_treat_dt(
    dt = df, 
    y_var = "first_Y", 
    group_var = "G", 
    id_var = "id", 
    birth_var = "birth_period", 
    hetero_var = "type"
)

summ_group_dt = create_group_first_treat_dt(
    summ_indiv_dt, 
    "first_Y", 
    "G", 
    c(1:time.periods), 
    unique(summ_indiv_dt$G),
    hetero_var = "type"
    )





quick_N_dt = create_N_per_period_from_summ(
    time_levels = unique(summ_group_dt$t), 
    summ_indiv_dt, 
    hetero_var = "type"
    )


test_that("Quick N per period works", {
    expect_equal(
        quick_N_dt,
        N_indiv_dt[, .(G, t, N, w, pr, type)]
    )
})



manual_infs = purrr::pmap(
    list(
        het_manual_fit$group, 
        het_manual_fit$time,
        het_manual_fit$type
    ),
    ~calculate_rc_influence_function(
        g_val = ..1, 
        t_val = ..2, 
        hetero_val = ..3,
        hetero_var = "type",
        lookup_indiv_table = summ_indiv_dt,
        row_id_var = "rowid",
        prop_score_known = TRUE
    )
)

clean_infs = function(inf_input) {
    new_infs = map(inf_input, "full_inf_func")
    new_adjustment = map(inf_input, "n_adjustment")
    new_infs = map2(new_infs, new_adjustment, ~.x*.y)
    inf_matrix = matrix(unlist(new_infs), nrow = nrow(new_infs[[1]]))
    return(inf_matrix)
}

manual_het_inf_matrix = clean_infs(manual_infs)


# Aggregate over type

combine_es_fun = function(et, att_df, inf_func) {
    w = att_df[, pr/sum(pr), .(event.time, treated)][event.time == et, V1]
    whichones = att_df[, event.time == et]
    weighted_het_if = inf_func[, whichones] %*% w
    return(weighted_het_if)
}

ets = het_manual_fit$event.time %>% unique() %>%
    sort()

combined_ifs = map(
    ets,
    ~combine_es_fun(.x, het_manual_fit, manual_het_inf_matrix)
)
et_se = map_dbl(
    combined_ifs,
    ~calculate_se(.x, biter = 1000, n_cores = 20, cluster_id = NULL)
)

unbalanced_panel_fit = did::att_gt(
    data = rc_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    est_method = "ipw",
    idname = "id",
    print_details = FALSE,
    panel = TRUE,
    allow_unbalanced_panel = TRUE,
    biters = 1000,
    control_group = "notyettreated" )

tidy_es_fit = unbalanced_panel_fit %>%
    aggte(type = "dynamic", biters = 1000) %>%
    tidy() %>%
    as_tibble()


comp_es_df = inner_join(
    tidy_es_fit,
    tibble(
        event.time = ets,
        manual_se = et_se
    ),
    by = "event.time"
)

test_that("Aggregated Het Ses roughly match unbalanced panel", {
    # Not the closest/most rigorous of checks
    diffs = comp_es_df %>% 
        filter(event.time >= 0) %>%
        mutate(diff = 100*abs(std.error - manual_se)/std.error) %>%
        pull()
    map(diffs, ~expect_lte(.x, 15))
})
