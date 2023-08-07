library(did)
library(BMisc)
library(dplyr)
library(purrr)
library(data.table)
library(furrr)
set.seed(10202)

ncl <- 1
time.periods <- 5
biters <- 200

custom_min = function(x) {if (length(x) > 0) min(x) else Inf}
# Creates simulation params
sim_params = did::reset.sim(
    time.periods = time.periods, 
    n = 2000)

sim_df = did::build_sim_dataset(sp_list = sim_params, panel = TRUE) %>%
    dplyr::as_tibble()

sim_df = sim_df %>%
    mutate(cluster = cluster %/% 3)

binary_sim_df = sim_df %>%
    group_by(G) %>%
    mutate(Y_above = Y > quantile(Y, 0.15)) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(
        first_Y =  custom_min(period[Y_above == TRUE]), 
        first_Y = if_else(!is.finite(first_Y), max(period), first_Y)
    ) %>%
    mutate(Y_binary = period >= first_Y)

df = as.data.table(binary_sim_df)

tictoc::tic()
rc_cs_fit = did::att_gt(
    data = binary_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    panel = FALSE,
    est_method = "ipw",
    idname = "id",
    print_details = FALSE,
    control_group = "notyettreated" )
tictoc::toc()



tictoc::tic()
cs_fit = did::att_gt(
    data = binary_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    est_method = "ipw",
    idname = "id",
    print_details = FALSE,
    bstrap = TRUE,
    biter = 10000,
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



#### Influence Function Time ####
N_indiv_dt = create_indiv_per_period_dt(df, "G", "period", c(1:time.periods), unique(df$G))
summ_indiv_dt = create_indiv_first_treat_dt(df, "first_Y", "G", "id")
summ_group_dt = create_group_first_treat_dt(
    summ_indiv_dt, 
    "first_Y", 
    "G", 
    c(1:time.periods), 
    unique(summ_indiv_dt$G))

manual_infs = purrr::map2(
    manual_did$group, 
    manual_did$time,
    ~calculate_influence_function(
        g_val = .x, 
        t_val = .y, 
        lookup_indiv_table = summ_indiv_dt,
        prop_score_known = FALSE
    )
)

new_infs = map(manual_infs, "full_inf_func")
new_adjustment = map(manual_infs, "n_adjustment")
new_infs = map2(new_infs, new_adjustment, ~.x*.y)
att_inf_output = cs_fit$inffunc
att_infs = map(1:ncol(att_inf_output), ~as.matrix(att_inf_output[, .x]))


inf_matrix = matrix(unlist(new_infs), nrow = nrow(new_infs[[1]]))
test_that("Influence functions match", {
    map2(new_infs, att_infs, all.equal) %>%
    map(expect_true)
}
)




test_that("Same analytical", {
    ed_V = Matrix::t(inf_matrix) %*% inf_matrix / nrow(inf_matrix)
    diff_df = bind_cols(
        ed_V = ed_V[, 1],
        cs_V = cs_fit$V_analytical[, 1]
    ) %>%
        mutate( 
            diff = ed_V - cs_V
        )
    map(diff_df$diff, expect_lte, 1e-6)
}
)



manual_se = calculate_se(
    inf_matrix,
    biter = 10000
    )

cs_manual_se = calculate_se(
    as.matrix(cs_fit$inffunc),
    biter = 10000
)



comp_se = tibble(
    manual_se = manual_se,
    cs_se = tidy_cs_fit %>%
        pull(std.error), 
    cs_manual_se = cs_manual_se
)


test_that("Homo SEs work", {
    # Results don't tend to match perfectly just due to bootstrapping it seems
    # Even taking cs_fit$inffunc and bootstrapping that tends to give slightly 
    # different results to using inf_matrix.
    se_diff = abs(tidy_cs_fit$std.error - manual_se)/tidy_cs_fit$std.error
    # no worse than 10% deviation from CS package
    map(100*se_diff, expect_lte, 10)
    # On average within 5%
    expect_lte(mean(se_diff), 0.05)

})


#### Aggregations ####

manual_es = manual_did %>%
    estimate_event_study(
        inf_matrix = inf_matrix, 
        y_var = "att_g_t",
        group_vector = summ_indiv_dt[, G],
        biter = 10000)
tidy_es_fit = cs_fit %>%
    aggte(type = "dynamic", biters = 10000) %>%
    tidy() %>%
    as_tibble()

tidy_es_fit %>%
    select(event.time, estimate, std.error)


test_that("Event Study Matches", {
    es_comp_estimate = bind_cols(
        cs_es = tidy_es_fit$estimate,
        manual_es = manual_es$estimate
    )
    es_comp_se = bind_cols(
        cs_es = tidy_es_fit$std.error,
        manual_es = manual_es$std.error
    )
    map2(
        es_comp_estimate$cs_es,
        es_comp_estimate$manual_es,
        ~expect_equal(.x, .y)
    )
    es_comp_se = es_comp_se %>%
        mutate(diff = abs(manual_es -cs_es)) %>%
        mutate(pct_diff = 100*diff/cs_es)
    map(
        es_comp_se$pct_diff,
        ~expect_lte(.x, 10)
    )

})

#### Group Averages ####
manual_group = manual_did %>%
    estimate_group_average(
        inf_matrix = inf_matrix, 
        y_var = "att_g_t",
        biter = 10000)
tidy_group_fit = cs_fit %>%
    aggte(type = "group", biters = 10000) %>%
    tidy() %>%
    as_tibble()

test_that("Group Aggregation Matches", {
    group_comp_estimate = bind_cols(
        cs_group = tidy_group_fit$estimate,
        manual_group = manual_group$estimate
    )
    group_comp_se = bind_cols(
        cs_group = tidy_group_fit$std.error,
        manual_group = manual_group$std.error
    )

    map2(
        group_comp_estimate$cs_group,
        group_comp_estimate$manual_group,
        ~expect_equal(.x, .y)
    )


    group_comp_se = group_comp_se %>%
        mutate(diff = abs(cs_group - manual_group)) %>%
        mutate(pct_diff = 100*diff/cs_group)
    map(
        group_comp_se$pct_diff,
        ~expect_lte(.x, 10)
    )

})

#### Clustered SEs ####

cluster_summ_indiv_dt = create_indiv_first_treat_dt(df, "first_Y", "G", "id", cluster_id = "cluster")

manual_cluster_se = calculate_se(
    inf_matrix,
    biter = 10000, 
    cluster_id = cluster_summ_indiv_dt[, cluster_id],
    cluster_id_2 = NULL
    )


tictoc::tic()
rc_cs_cluster_fit = did::att_gt(
    data = binary_sim_df,
    yname = "Y_binary",
    tname = "period",
    gname = "G",
    panel = FALSE,
    clustervars = "cluster",
    est_method = "ipw",
    idname = "id",
    print_details = FALSE,
    control_group = "notyettreated" )
tictoc::toc()




test_that("Clustered SEs work", {
    comp_se_df = bind_cols(
        manual = manual_cluster_se, 
        cs = rc_cs_cluster_fit %>%
            tidy() %>%
            pull(std.error)
    ) %>%
        mutate(diff = 100*(cs - manual)/cs) 

    # Not more than 10% smaller off for any indiv error
    map(
        comp_se_df$diff, 
        expect_lte, 
        10
    )    

    # on average not off by more than 5%
    expect_lte(mean(comp_se_df$diff), 5)

})

