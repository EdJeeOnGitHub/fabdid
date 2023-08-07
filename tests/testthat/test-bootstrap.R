library(BMisc)
library(data.table)
library(dplyr)
library(purrr)

set.seed(1232)
N_indiv = 52
N_att = 200
if_test = matrix(rnorm(N_indiv*N_att), nrow = N_indiv, ncol = N_att)


test_that("IF output exists", {
    bc_if = run_multiplier_bootstrap(if_test, 100)
    ed_if = run_nested_multiplier_bootstrap(if_test, 100)

    expect_equal(dim(ed_if), dim(bc_if))
})


test_that("R Faster", {
    mult_boot_bench = microbenchmark::microbenchmark(
        bc = run_multiplier_bootstrap(if_test, 100),
        ed = run_nested_multiplier_bootstrap(if_test, 100)
    )
    bc_time = mean(mult_boot_bench$time[mult_boot_bench$expr == "bc"])
    ed_time = mean(mult_boot_bench$time[mult_boot_bench$expr == "ed"])
    expect_lte(ed_time, bc_time)
})


test_that("Standard Errors Match-ish", {
    bc_se = calculate_se(if_test, bs_fun_type = "bc")
    ed_se = calculate_se(if_test, bs_fun_type = "ed")
    se_diff = abs(bc_se - ed_se)
    expect_true(all(se_diff < 1e-2))
})



test_that("R SEs Faster", {
    se_bench = microbenchmark::microbenchmark(
        bc_se = calculate_se(if_test, bs_fun_type = "bc"),
        ed_se = calculate_se(if_test, bs_fun_type = "ed"),
        times = 10
    )
    bc_time = mean(se_bench$time[se_bench$expr == "bc_se"])
    ed_time = mean(se_bench$time[se_bench$expr == "ed_se"])
    expect_lte(ed_time, bc_time)
})


## Now test clustered standard errors
N_clusters = 12
cluster_id = sample(1:N_clusters, replace = TRUE, N_indiv)

clust_if = run_nested_multiplier_bootstrap(if_test, 1, cluster_id_2 = cluster_id)
clust_se = calculate_se(
    if_test,
    bs_fun_type = "ed",
    cluster_id_2 = cluster_id
)

test_that("Clustered SEs Exist",  {
    expect_equal(length(clust_se), N_att)
}
)
