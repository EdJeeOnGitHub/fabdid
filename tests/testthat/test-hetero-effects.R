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
    n = 10000)

probs = c(1:10)/sum(1:10)

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
        type = sample(1:10, size = 1, prob = probs)
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


manual_homo_inf_matrix = clean_infs(no_het_fit$inf_func_output)

## Het ES
het_manual_es = het_manual_did %>%
    estimate_event_study(
        inf_matrix = manual_het_inf_matrix, 
        y_var = "att_g_t",
        group_vector = summ_indiv_dt[, G],
        prop_score_known = TRUE,
        biter = 10000, 
        hetero_var = "type", 
        n_cores = 2
        )

# Testing this vs package kinda hard. Need to write DGP to test this

agg_het_es = het_manual_did %>%
    estimate_event_study(
        inf_matrix = manual_het_inf_matrix, 
        y_var = "att_g_t",
        # group_vector = summ_indiv_dt[, G],
        prop_score_known = TRUE,
        biter = 1000, 
        n_cores = 2
        )


het_es = het_manual_did %>%
    estimate_event_study(
        inf_matrix = manual_het_inf_matrix, 
        y_var = "att_g_t",
        hetero_var = "type",
        # group_vector = summ_indiv_dt[, G],
        prop_score_known = TRUE,
        biter = 1000, 
        n_cores = 2
        )

test_that(
    "IF standard errors for hetero effects match theory heuristic", {

        comp_es_size = copy(het_es)
        comp_es_size = merge(
            comp_es_size,
            agg_het_es[, .(full_std_error = std.error, event.time)],
            by  = "event.time"
        )

        prob_df = tibble(
            type = 1:10,
            probs = probs
        )

        comp_es_size = merge(
            comp_es_size,
            prob_df,
            by = "type"
        )
        comp_es_size[, emp_ratio := std.error/full_std_error]
        comp_es_size[, theory_ratio := 1/sqrt(probs)]

        # this shouldn't be perfect since this isn't just comparing std.error of means with 
        # smaller samples (for example, nyt always present + IF accounts for correlation across 
        # thetas) but as a rule of thumb we hope this is between 0.75-1.25
        # comp_es_size %>%
        #     filter(event.time > 0) %>%
        #     ggplot(aes( 
        #         x = theory_ratio,
        #         y = emp_ratio
        #     )) +
        #     geom_point() +
        #     geom_abline()

        emp_theory_relationship = comp_es_size %>%
            filter(event.time > 0) %>%
            lm(
                emp_ratio ~  theory_ratio,
                data = .
            ) %>%
            tidy(conf.int = TRUE) %>%
            filter(term == "theory_ratio")  %>%
            pull(estimate)


        expect_lte(emp_theory_relationship, 1.25)
        expect_gte(emp_theory_relationship, 0.75)
    }
)


combine_es_fun = function(et, att_df, inf_func) {
    w = att_df[, pr/sum(pr), .(event.time, treated)][event.time == et, V1]
    whichones = att_df[, event.time == et]
    weighted_het_if = inf_func[, whichones] %*% w
    return(weighted_het_if)
}


ets = het_manual_did$event.time %>% unique() %>% sort()
combined_ifs = map(
    ets,
    ~combine_es_fun(.x, het_manual_did, manual_het_inf_matrix)
)

weighted_manual_homo_inf_matrix = do.call(cbind, combined_ifs)


et_se = map_dbl(
    combined_ifs,
    ~calculate_se(.x, biter = 1000, n_cores = 2, cluster_id = NULL)
)

test_that(
    "Combining IFs recovers event study SEs", {
        se_diff = 100*(et_se - agg_het_es$std.error)/agg_het_es$std.error
        map(se_diff, ~expect_lte(abs(.x), 15))
    }
)




agg_manual_es = no_het_fit$att_df %>%
    estimate_event_study(
        inf_matrix = manual_homo_inf_matrix, 
        y_var = "att_g_t",
        prop_score_known = TRUE,
        biter = 10000, 
        n_cores = 2
        )


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








