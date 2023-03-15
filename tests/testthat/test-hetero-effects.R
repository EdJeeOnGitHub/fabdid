
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
        first_Y =  min(period[Y_above == TRUE]), 
        first_Y = if_else(!is.finite(first_Y), max(period), first_Y)
    ) %>%
    mutate(Y_binary = period >= first_Y) %>%
    group_by(id) %>%
    mutate(
        type = sample(1:2, size = 1)
    )


df = as.data.table(binary_sim_df)


df[, .N, .(G)]

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

no_het_fit = estimate_did(
    df,
    y_var = "first_Y",
    group_var = "G",
    t_var = "period",
    id_var = "id"
)$att_df


het_manual_did = estimate_did(
    df,
    y_var = "first_Y",
    group_var = "G",
    t_var = "period",
    id_var = "id",
    hetero_var = "type"
)$att_df

manual_did = het_manual_did[, .(att_g_t = sum(att_g_t*pr_hetero)), .(group, time)]

comp_df = inner_join(
    manual_did %>% rename(manual_estimate = att_g_t),
    tidy_cs_fit %>%
        select(group, time, cs_estimate = estimate), 
    by = c("group","time")
) %>%
    left_join(
        no_het_fit %>% rename(no_het_manual_estimate = att_g_t),
        by = c("group", "time")
    )

no_het_lookup_max_error = comp_df %>%
    mutate(error = abs(no_het_manual_estimate - cs_estimate)) %>%
    summarise(max_error = max(error)) %>%
    pull()
lookup_max_error = comp_df %>%
    mutate(error = abs(manual_estimate - cs_estimate)) %>%
    summarise(max_error = max(error)) %>%
    pull()

test_that("Manual Estimates OK", {
    expect_lte(lookup_max_error, 1e-8)
    expect_lte(no_het_lookup_max_error, 1e-8)
}
)

