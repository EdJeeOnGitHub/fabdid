
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
    n = 50)

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
    mutate(Y_binary = period >= first_Y) %>%
    group_by(id) %>%
    mutate(
        type = sample(1:2, size = 1, prob = c(0.7, 0.3))
    )



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

no_het_fit = estimate_did(
    df,
    y_var = "first_Y",
    group_var = "G",
    t_var = "period",
    id_var = "id",
    prop_score_known = TRUE
)

no_het_manual_did = no_het_fit$att_df

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
        no_het_manual_did %>% rename(no_het_manual_estimate = att_g_t),
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



#### Influence Function Time ####
N_indiv_dt = create_indiv_per_period_dt(
    df, 
    "G", 
    "period", 
    c(1:time.periods), 
    unique(df$G), 
    hetero_var = "type")
summ_indiv_dt = create_indiv_first_treat_dt(
    df, 
    "first_Y", 
    "G", 
    "id", 
    hetero_var = "type")
summ_group_dt = create_group_first_treat_dt(
    summ_indiv_dt, 
    "first_Y", 
    "G", 
    c(1:time.periods), 
    unique(summ_indiv_dt$G), 
    hetero_var = "type"
    )


manual_infs = purrr::pmap(
    list(
        het_manual_did$group, 
        het_manual_did$time,
        het_manual_did$type
    ),
    ~calculate_influence_function(
        g_val = ..1, 
        t_val = ..2, 
        hetero_val = ..3,
        hetero_var = "type",
        lookup_indiv_table = summ_indiv_dt,
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
att_inf_output = cs_fit$inffunc
att_infs = map(1:ncol(att_inf_output), ~as.matrix(att_inf_output[, .x]))

test_that("Influence function has same dim as att_df", {
    expect_equal(ncol(manual_het_inf_matrix), nrow(het_manual_did))
}
)
#### Testing we can recover IF

manual_infs = purrr::pmap(
    list(
        het_manual_did$group, 
        het_manual_did$time,
        het_manual_did$type
    ),
    ~calculate_influence_function(
        g_val = ..1, 
        t_val = ..2, 
        hetero_val = ..3,
        hetero_var = "type",
        lookup_indiv_table = summ_indiv_dt,
        prop_score_known = TRUE
    )
)


homo_2_2_if = calculate_influence_function(
    g_val = 2,
    t_val = 2,
    hetero_val = NULL,
    hetero_var = NULL,
    lookup_indiv_table = summ_indiv_dt,
    prop_score_known = TRUE
)


hetero_2_2_1_if = calculate_influence_function(
    g_val = 2,
    t_val = 2,
    hetero_val = 1,
    hetero_var = "type",
    lookup_indiv_table = summ_indiv_dt,
    prop_score_known = TRUE
)

hetero_2_2_2_if = calculate_influence_function(
    g_val = 2,
    t_val = 2,
    hetero_val = 2,
    hetero_var = "type",
    lookup_indiv_table = summ_indiv_dt,
    prop_score_known = TRUE
)

hetero_2_2_1_if
hetero_2_2_2_if


a = 0.2
b = 1 - a
z = rnorm(1000)
x = z + rnorm(1000)
y = a*x + b*z
if_z = z - mean(z)
if_x = x - mean(x)

if_y = y - mean(y)
if_comp = a*if_x + b*if_z

if_y - if_comp


het_weight_dt = het_manual_did[, .(type, weight = pr/sum(pr)), .(group, time)] %>%
    tidyr::spread(type, weight, sep = "_")

create_homo_if = function(g_val,
                          t_val,
                          lookup_indiv_table, 
                          het_att_df) {
    if_1 = calculate_influence_function(
        g_val = g_val,
        t_val = t_val,
        hetero_val = 1,
        hetero_var = "type",
        lookup_indiv_table = lookup_indiv_table,
        prop_score_known = TRUE
    )
    if_2 = calculate_influence_function(
        g_val = g_val,
        t_val = t_val,
        hetero_val = 2,
        hetero_var = "type",
        lookup_indiv_table = lookup_indiv_table,
        prop_score_known = TRUE
    )
    het_weight_dt = het_att_df[, .(type, weight = pr/sum(pr)), .(group, time)] %>%
        tidyr::spread(type, weight, sep = "_")
    w_1 = het_weight_dt[group == g_val & time == t_val, type_1]
    w_2 = het_weight_dt[group == g_val & time == t_val, type_2]
    combined_if = w_1*if_1$full_inf_func*if_1$n_adjustment + w_2*if_2$full_inf_func*if_2$n_adjustment
    homo_if = calculate_influence_function(
        g_val = g_val,
        t_val = t_val,
        lookup_indiv_table = lookup_indiv_table,
        prop_score_known = TRUE
    )
    return(lst(
        if_1,
        if_2,
        combined_if,
        homo_if
    ))
}


comp_if_2_2 = create_homo_if(2, 2, summ_indiv_dt, het_manual_did)

comp_df = tibble(
    homo = comp_if_2_2$homo_if$full_inf_func[, 1],
    combined = comp_if_2_2$combined_if[, 1]
) %>%
    mutate(
        rowid = 1:n()
    )

comp_if_2_2$if_1$n_adjustment*comp_if_2_2$if_2$n_adjustment

comp_df %>%
    mutate(
        pct_diff = homo/combined
    )  %>%
    pull(pct_diff) %>%
    unique()

comp_df

problem_rowids = comp_df %>%
    filter(abs(homo - combined) > 1e-3) %>%
    pull(rowid)
problem_rowids
comp_df %>%
    filter(homo != combined)
comp_if_2_2$if_2$full_inf_func[problem_rowids, ]

summ_indiv_dt %>%
    filter(rowid %in% problem_rowids) %>%
    select(-const) %>%
    select(-id)

summ_indiv_dt %>%
    filter(!(rowid %in% problem_rowids)) %>%
    select(-const) %>%
    select(-id)

comp_df %>%
    filter(rowid %in% problem_rowids)


stop()
test_calculate_influence_function(
        g_val = 2,
        t_val = 2,
        hetero_val = 2,
        hetero_var = "type",
        lookup_indiv_table = summ_indiv_dt,
        prop_score_known = TRUE
    )


test_calculate_influence_function(
        g_val = 2,
        t_val = 2,
        lookup_indiv_table = summ_indiv_dt,
        prop_score_known = TRUE
    )

comp_df %>%
    mutate(
        rowid = 1
    )

summ_indiv_dt[G == 2]

comp_df %>%
    mutate(diff = homo - combined)

ed = cbind(
)

ed

comp_if_2_2$if_1$n_adjustment
1 - comp_if_2_2$if_2$n_adjustment

(ed[, 1]/ed[, 2]) %>% unique()
(ed[, 2]/ed[, 1]) %>% unique()


manual_homo_inf_matrix = clean_infs(no_het_fit$inf_func_output)
combine_fun = function(group_val, time_val, att_df, inf_func) {
    w = att_df[, pr/sum(pr), .(group, time)][group == group_val & time == time_val, V1]
    whichones = att_df[, group == group_val & time == time_val]
    weighted_if = inf_func[, whichones] %*% w
    return(weighted_if)
}

comp_fun = function(group_val, time_val, homo_att_df, homo_inf_func, het_att_df, het_inf_func) {
    combined_if = combine_fun(group_val, time_val, het_att_df, het_inf_func)
    which_idx = homo_att_df[, group == group_val & time == time_val]
    homo_if = homo_inf_func[, which_idx]
    return(lst(homo_if, combined_if[, 1]))
}



ed = comp_fun(
    2,
    2,
    manual_did,
    manual_homo_inf_matrix,
    het_manual_did,
    manual_het_inf_matrix
)

ed

df = tibble(
    homo = ed$homo_if,
    combined = ed$combined_if
)


df

df %>%
    summarise_all(var)

df %>%
    ggplot(aes(
        x = homo, 
        y = combined
    )) +
    geom_point() +
    geom_abline()

manual_het_inf_matrix[3, ]
manual_homo_inf_matrix[3, ]

binary_sim_df

df





## Het ES
het_manual_es = het_manual_did %>%
    estimate_event_study(
        inf_matrix = manual_het_inf_matrix, 
        y_var = "att_g_t",
        group_vector = summ_indiv_dt[, G],
        prop_score_known = TRUE,
        biter = 10000, 
        hetero_var = "type")

# Testing this vs package kinda hard. Need to write DGP to test this
test_that("Correct number of rows", {
    expect_equal(nrow(het_manual_es), 14)
})

agg_het_es = het_manual_did %>%
    estimate_event_study(
        inf_matrix = manual_het_inf_matrix, 
        y_var = "att_g_t",
        # group_vector = summ_indiv_dt[, G],
        prop_score_known = TRUE,
        biter = 1000)


combine_fun(2, 2, het_manual_did, manual_het_inf_matrix)





combine_es_fun = function(et, att_df, inf_func) {
    w = att_df[, pr/sum(pr), .(event.time, treated)][event.time == et, V1]
    whichones = att_df[, event.time == et]
    weighted_het_if = inf_func[, whichones] %*% w
    return(weighted_het_if)
}


ets = het_manual_did$event.time %>% unique() %>% sort()
combined_ifs = map(
    ets,
    ~combine_fun(.x, het_manual_did, manual_het_inf_matrix)
)

weighted_manual_homo_inf_matrix = do.call(cbind, combined_ifs)
manual_homo_inf_matrix


et_se = map_dbl(
    combined_ifs,
    ~calculate_se(.x, biter = 1000, n_cores = 20, cluster_id = NULL)
)

et_se

agg_het_es$ed = et_se

agg_het_es[, .(std.error, ed)]



str(weighted_manual_homo_inf_matrix)


ed = combine_fun(0, het_manual_did, manual_het_inf_matrix)

weighted_het_if = manual_het_inf_matrix %*% w_het

str(weighted_het_if)
str(manual_het_inf_matrix)


manual_es = no_het_fit$att_df %>%
    estimate_event_study(
        inf_matrix = manual_homo_inf_matrix, 
        y_var = "att_g_t",
        group_vector = summ_indiv_dt[, G],
        biter = 10000)


agg_manual_es = no_het_fit$att_df %>%
    estimate_event_study(
        inf_matrix = manual_homo_inf_matrix, 
        y_var = "att_g_t",
        prop_score_known = TRUE,
        biter = 10000)

agg_het_es
agg_manual_es

## Test that aggregating over het effects event study gives similar results to 
# running unconditional event study


tidy_es_fit = cs_fit %>%
    aggte(type = "dynamic", biters = 10000) %>%
    tidy() %>%
    as_tibble()

test_that("Het and Agg Event Study Matches", {
    es_comp_estimate = bind_cols(
        agg_het_es = agg_het_es$estimate,
        cs_es = tidy_es_fit$estimate
    )
    es_comp_se = bind_cols(
        agg_het_es = agg_het_es$std.error,
        cs_es = tidy_es_fit$std.error
    )
    map2(
        es_comp_estimate$agg_het_es,
        es_comp_estimate$cs_es,
        ~expect_equal(.x, .y)
    )
    es_comp_se = es_comp_se %>%
        mutate(diff = abs(cs_es - agg_het_es)) %>%
        mutate(pct_diff = 100*diff/cs_es)
    map(
        es_comp_se$pct_diff,
        ~expect_lte(.x, 10)
    )

})


agg_manual_es %>%
    select(event.time, estimate, std.error) %>%
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

