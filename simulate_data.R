# Intra-articular corticosteroid injections for symptomatic osteoarthritis:
# The effects of injection dose, frequency, and timing on safety
#
# Simulate data with which to develop aim 3 analysis code
#
# References
# [1] [Generating survival times to simulate Cox proportional hazards models](https://doi.org/10.1002/sim.2059)
# [2] [Simulating survival outcomes: setting the parameters for the desired distribution](https://www.rdatagen.net/post/2022-02-08-simulating-survival-outcomes-setting-the-parameters-for-the-desired-distribution/)
# [3] [Follow-up: simstudy function for generating parameters for survival distribution](https://www.rdatagen.net/post/2022-02-22-follow-up-simstudy-function-for-generating-parameters-for-survival-distribution/)
# [4] [Gemini query:
# How can I generate simulated survival times from a distribution defined by survival at two specific points during follow up? For example, I might know that survival at time 100 is 80% and that survival at time 200 is 20%. I want to code this in R.
# ]
#
# Date created: 20261111
# Last modified: 20261115
# John M Sahrmann


# * Setup

library(dplyr)
library(purrr)


# * Constant definitions


# * Function definitions


# ** Person-level generators, i.e., generate one value for each of `n` patients

# Generate event or censoring times from a discrete distribution with
# specific event probabilities at particular event times.
gen_discrete_times <- function(n, event_probs, event_times) {
  if (!any(is.na(event_times))) {
    event_times <- c(event_times, NA)
    event_probs <- c(event_probs, 1 - sum(event_probs))
  }
  sample(event_times, size = n, replace = TRUE, prob = event_probs)
}

# Generate event or censoring times from an exponential distribution
# with survival S1 at time t1.  Inspired by [2] and [3] with guidance
# from Gemini [4].
gen_exp_times <- function(n, S1, t1) {
  rate <- -log(S1) / t1
  rexp(n, rate = rate)
}

# Generate event or censoring times from a uniform distribution with
# survival S1 at time t1.
gen_unif_times <- function(n, S1, t1) {
  max <- t1 / (1 - S1)
  runif(n, min = 0, max = max)
}

# Generate event or censoring times from a Weibull distribution with
# survival S1, S2 at times t1, t2.  Inspired by [2] and [3] with
# guidance from Gemini [4].
gen_weib2_times <- function(n, S1, t1, S2, t2) {
  k <- log(log(S2) / log(S1)) / log(t2 / t1)
  lambda <- t1 / (-log(S1))^(1 / k)
  u <- runif(n)
  lambda * (-log(u))^(1 / k)
}

# Generate event or censoring times from a Weibull distribution with
# (approximate) survival at times specified in a list of the form
# ((S1, t1), (S2, t2), ..., (Sn, tn)).  Adapted from [2] and [3].
gen_weibn_times <- function(n, points) {
  loss_weib <- function(params, points) {
    loss_weib_1 <- function(point) {
      (params[[2]] * (log(-log(point[[1]])) - params[[1]]) - log(point[[2]]))^2
    }
    sum(sapply(points, \(p) loss_weib_1(p)))
  }
  r <- optim(
    par = c(1, 1), fn = loss_weib, points = points,
    method = "L-BFGS-B", lower = c(-Inf, 0), upper = c(Inf, Inf)
  )
  f <- r$par[[1]]; shape <- r$par[[2]]
  u <- runif(n)
  (-(log(u) / exp(f)))^shape
}

# Generate ages from a normal distribution with a given mean and
# standard deviation, but keeping only ages greater than or equal to
# some minimum.
gen_ages <- function(n, mean, sd, min) {
  ages <- numeric(n)
  i <- 1
  while (i <= n) {
    age <- rnorm(1, mean, sd)
    if (age >= min) {
      ages[[i]] <- age
      i <- i + 1
    }
  }
  ages
}


# ** Person-time-level generators, i.e., generate `end_fu` values for each of `n` patients

# Define event indicators for each time point from 1 to `end_fu` for
# `n` event times drawn from an exponential distribution with survival
# S1 at time t1.  Before events, the values are zero, at the times of
# events, the values are one, and after events, the values are `NA`.
gen_exp_series <- function(n, end_fu, S1, t1) {
  times <- ceiling(gen_exp_times(n = n, S1 = S1, t1 = t1))
  are_events <- times <= end_fu
  series_matrix <- matrix(nrow = end_fu, ncol = n)
  for (i in seq_len(n)) {
    if (are_events[[i]]) {
      series <- c(
        rep(0, times = times[[i]] - 1),
        1,
        rep(NA, times = end_fu - times[[i]])
      )
    } else {
      series <- rep(0, times = end_fu)
    }
    series_matrix[, i] <- series
  }
  as.integer(series_matrix)
}

# Define event indicators for each time point from 1 to `end_fu` for
# `n` event times drawn from a uniform distribution with survival S1
# at time t1.  Before events, the values are zero, at the times of
# events, the values are one, and after events, the values are `NA`.
gen_unif_series <- function(n, end_fu, S1, t1) {
  times <- ceiling(gen_unif_times(n = n, S1 = S1, t1 = t1))
  are_events <- times <= end_fu
  series_matrix <- matrix(nrow = end_fu, ncol = n)
  for (i in seq_len(n)) {
    if (are_events[[i]]) {
      series <- c(
        rep(0, times = times[[i]] - 1),
        1,
        rep(NA, times = end_fu - times[[i]])
      )
    } else {
      series <- rep(0, times = end_fu)
    }
    series_matrix[, i] <- series
  }
  as.numeric(series_matrix)
}

# Define event indicators for each time point from 1 to `end_fu` for
# `n` event times drawn from a Weibull distribution with survival S1,
# S2 at times t1, t2.  Before events, the values are zero, at the
# times of events, the values are one, and after events, the values
# are `NA`.
gen_weib2_series <- function(n, end_fu, S1, t1, S2, t2) {
  times <- ceiling(gen_weib2_times(n = n, S1 = S1, t1 = t1, S2 = S2, t2 = t2))
  are_events <- times <= end_fu
  series_matrix <- matrix(nrow = end_fu, ncol = n)
  for (i in seq_len(n)) {
    if (are_events[[i]]) {
      series <- c(
        rep(0, times = times[[i]] - 1),
        1,
        rep(NA, times = end_fu - times[[i]])
      )
    } else {
      series <- rep(0, times = end_fu)
    }
    series_matrix[, i] <- series
  }
  as.integer(series_matrix)
}

# Define frailty series from 1 to `end_fu` for `n` patients.  Patients
# are classified with equal probability into one of two
# phenotypes---low frailty with initial frailty drawn from N(0.17,
# 0.01) and shocks generated from N(0, 0.01), and high frailty with
# initial frailty drawn from N(0.30, 0.02) and shocks generated from
# N(0, 0.02).  Values after the initial point are drawn from an AR(1)
# model with parameter 0.6.
gen_frailty_series <- function(n, end_fu) {
  phi <- 0.6
  series_matrix <- matrix(nrow = end_fu, ncol = n)
  for (i in seq_len(n)) {
    series <- numeric(end_fu)
    frailty_type <- sample(c("low frailty", "high frailty"), size = 1, prob = c(0.5, 0.5))
    if (frailty_type == "low frailty") {
      series[[1]] <- rnorm(1, 0.17, 0.01)
      for (t in 2:end_fu) {
        series[[t]] <- series[[1]] + phi * (series[[t-1]] - series[[1]]) + rnorm(1, 0, 0.01)
      }
      series_matrix[, i] <- series
    }
    else if (frailty_type == "high frailty") {
      series[[1]] <- rnorm(1, 0.30, 0.02)
      for (t in 2:end_fu) {
        series[[t]] <- series[[1]] + phi * (series[[t-1]] - series[[1]]) + rnorm(1, 0, 0.02)
      }
      series_matrix[, i] <- series
    }
  }
  as.numeric(series_matrix)
}

# Generate binary series from 1 to `end_fu` for `n` patients.  The
# number of 'events' for each patient is drawn from Poisson(`mean`),
# then event times are drawn without replacement from 1 to `end_fu`.
# Series values are one at event times and the two time points after
# events, and zero otherwise.
gen_time_varying_series <- function(n, end_fu, mean) {
  series_matrix <- matrix(nrow = end_fu, ncol = n)
  for (i in seq_len(n)) {
    series <- numeric(end_fu)
    n_events <- rpois(1, mean)
    event_times <- sample(seq_len(end_fu), size = n_events)
    series[event_times] <- 1
    series <- pmax(series, lag(series, 1), lag(series, 2), na.rm = TRUE)
    series_matrix[, i] <- series
  }
  as.integer(series_matrix)
}

# Generate binary series from 1 to `end_fu` for `n` patients.  A
# single 'event' time is drawn for each patient.  Series values are
# one at event times and all time points thereafter, and zero
# otherwise.
gen_step_series <- function(n, end_fu, S1, t1, S2, t2) {
  times <- ceiling(gen_weib2_times(n = n, S1 = S1, t1 = t1, S2 = S2, t2 = t2))
  are_events <- times <= end_fu
  series_matrix <- matrix(nrow = end_fu, ncol = n)
  for (i in seq_len(n)) {
    if (are_events[[i]]) {
      series <- c(
        rep(0, times = times[[i]] - 1),
        1,
        rep(1, times = end_fu - times[[i]])
      )
    } else {
      series <- rep(0, times = end_fu)
    }
    series_matrix[, i] <- series
  }
  as.integer(series_matrix)
}


# ** Data set generator

#' Simulate data for aim 3 analyses
#'
#' We generate data for `n` patients using named lists of
#' data-generating functions.  Each function in
#' `time_fixed_generator_list` should be a function of one variable
#' that produces `n` observations.  Each function in
#' `time_varying_generator_list` and `time_to_event_generator_list`
#' should be a function of two variables that produces `end_fu`
#' observations for each patient (for `n` * `end_fu` observations
#' total).  (`time_varying_generator_list` and
#' `time_to_event_generator_list` act identically but are separated
#' for conceptual clarity.)
#'
#' Once the data are generated, they're merged to form a single data
#' set.  If only time-fixed generators are given, the output data set
#' will contain `n` observations; if either of
#' `time_varying_generator_list` or `time_to_event_generator_list` are
#' given, the output data set will contain `n` * `end_fu`
#' observations.  An additional column `id` with values from 1 to `n`
#' is automatically added as the first column.  If time-varying or
#' time-to-event data are created, a column `time` with values from 1
#' to `end_fu` repeated for each patient is added before the
#' time-varying / time-to-event data columns.
#'
#' @param n The number of patients.
#' @param time_fixed_generator_list A named list of functions to
#'   generate patient-level data.
#' @param time_varying_generator_list A named list of functions to
#'   generate patient-time-level data.
#' @param time_to_event_generator_list A named list of functions to
#'   generate patient-time-level data.
#' @param end_fu The length of follow up if time-varying or
#'   time-to-event data are being generated.  Time is measured from
#'   one to `end_fu`.
#'
#' @return A Tibble with columns `id`, `time` (if either
#'   `time_varying_generator_list` or `time_to_event_generator_list`
#'   are not `NULL`), and one column of data for each function in the
#'   data generator lists with names matching those in the
#'   corresponding list.
simulate_data <- function(
  n,
  time_fixed_generator_list = NULL,
  time_varying_generator_list = NULL,
  time_to_event_generator_list = NULL,
  end_fu = NULL
) {
  if (!is.null(time_fixed_generator_list)) {
    time_fixed_data <- purrr::map(time_fixed_generator_list, \(fn) fn(n)) |>
      tibble::as_tibble() |>
      dplyr::mutate(id = 1:n, .before = 1)
  }
  if (!is.null(time_varying_generator_list)) {
    time_varying_data <- purrr::map(time_varying_generator_list, \(fn) fn(n, end_fu)) |>
      tibble::as_tibble() |>
      dplyr::mutate(
        id = rep(1:n, each = end_fu),
        time = rep(1:end_fu, times = n),
        .before = 1
      )
  }
  if (!is.null(time_to_event_generator_list)) {
    time_to_event_data <- purrr::map(time_to_event_generator_list, \(fn) fn(n, end_fu)) |>
      tibble::as_tibble() |>
      dplyr::mutate(
        id = rep(1:n, each = end_fu),
        time = rep(1:end_fu, times = n),
        .before = 1
      )
  }
  if (is.null(time_fixed_generator_list) && is.null(time_varying_generator_list) && !is.null(time_to_event_generator_list)) {
    merged_data <- time_to_event_data
  } else if (is.null(time_fixed_generator_list) && !is.null(time_varying_generator_list) && is.null(time_to_event_generator_list)) {
    merged_data <- time_varying_data
  } else if (is.null(time_fixed_generator_list) && !is.null(time_varying_generator_list) && !is.null(time_to_event_generator_list)) {
    merged_data <- dplyr::inner_join(time_varying_data, time_to_event_data, by = dplyr::join_by(id, time))
  } else if (!is.null(time_fixed_generator_list) && is.null(time_varying_generator_list) && is.null(time_to_event_generator_list)) {
    merged_data <- time_fixed_data
  } else if (!is.null(time_fixed_generator_list) && is.null(time_varying_generator_list) && !is.null(time_to_event_generator_list)) {
    merged_data <- dplyr::inner_join(time_fixed_data, time_to_event_data, by = dplyr::join_by(id))
  } else if (!is.null(time_fixed_generator_list) && !is.null(time_varying_generator_list) && is.null(time_to_event_generator_list)) {
    merged_data <- dplyr::inner_join(time_fixed_data, time_varying_data, by = dplyr::join_by(id))
  } else if (!is.null(time_fixed_generator_list) && !is.null(time_varying_generator_list) && !is.null(time_to_event_generator_list)) {
    merged_data <- dplyr::inner_join(time_fixed_data, time_varying_data, by = dplyr::join_by(id)) |>
      dplyr::inner_join(time_to_event_data, by = dplyr::join_by(id, time))
  }
  merged_data
}


# * Scratch work


# ** Example 1, time-fixed-only with base R functions

set.seed(56485463)

ex1_fixed_generators <- list(
  y = \(n) rnorm(n),
  x = \(n) rbinom(n, size = 1, prob = 0.5)
)

ex1_data <- simulate_data(
  n = 20,
  time_fixed_generator_list = ex1_fixed_generators
)
ex1_data
#> # A tibble: 20 × 3
#>       id       y     x
#>    <int>   <dbl> <int>
#>  1     1 -1.87       0
#>  2     2 -0.890      0
#>  3     3  0.809      1
#>  4     4 -1.21       1
#>  5     5  1.66       0


# ** Example 2, time-fixed-only with custom functions

set.seed(13495295)

ex2_fixed_generators <- list(
  age = \(n) gen_ages(n = n, mean = 37, sd = 13, min = 67),
  sex = \(n) sample(c("F", "M"), size = n, replace = TRUE, prob = c(0.5, 0.5)),
  # 30% mortality at 60 months = 5 years
  death_time = \(n) gen_exp_times(n = n, S1 = 0.7, t1 = 60),
  # 40% loss to follow up at 5 years
  censor_time = \(n) gen_unif_times(n = n, S1 = 0.6, t1 = 60)
)

ex2_data <- simulate_data(
  n = 100,
  time_fixed_generator_list = ex2_fixed_generators
)
ex2_data
#> # A tibble: 100 × 5
#>       id   age sex   death_time censor_time
#>    <int> <dbl> <chr>      <dbl>       <dbl>
#>  1     1  77.8 F          138.       24.0  
#>  2     2  76.7 M           15.8       0.583
#>  3     3  72.5 F          261.       45.5  
#>  4     4  70.1 M           38.4      74.4  
#>  5     5  67.7 F           31.2      93.0  

mean(ex2_data$death_time < 60)
#> [1] 0.3
mean(ex2_data$censor_time < 60)
#> [1] 0.37


# ** Example 3, time-varying-only

set.seed(4662765)

ex3_varying_generators <- list(
  frailty = \(n, end_fu) gen_frailty_series(n = n, end_fu = end_fu)
)

ex3_data <- simulate_data(
  n = 10,
  time_varying_generator_list = ex3_varying_generators,
  end_fu = 20
)
ex3_data
#> # A tibble: 200 × 3
#>       id  time frailty
#>    <int> <int>   <dbl>
#>  1     1     1   0.178
#>  2     1     2   0.186
#>  3     1     3   0.181
#>  4     1     4   0.177
#>  5     1     5   0.176
#>         [...]
#> 16     1    16   0.159
#> 17     1    17   0.165
#> 18     1    18   0.170
#> 19     1    19   0.161
#> 20     1    20   0.178
#> 21     2     1   0.275
#> 22     2     2   0.282
#> 23     2     3   0.283
#> 24     2     4   0.291
#> 25     2     5   0.277


# ** Example 4, time-fixed, time-varying, and time-to-event

set.seed(90891157)

ex4_fixed_list <- list(
  g = \(n) rpois(n, lambda = 0.5),
  x = \(n) rbinom(n, size = 1, prob = 0.5)
)
ex4_varying_list <- list(
  f = \(n, end_fu) gen_frailty_series(n = n, end_fu = end_fu),
  l = \(n, end_fu) gen_step_series(n = n, end_fu = end_fu, S1 = 0.5, t1 = 1, S2 = 0.3, t2 = 10)
)
ex4_event_list <- list(
  y = \(n, end_fu) gen_exp_series(n = n, end_fu = end_fu, S1 = 0.75, t1 = 10)
)

ex4_data <- simulate_data(
  n = 20,
  time_fixed_generator_list = ex4_fixed_list,
  time_varying_generator_list = ex4_varying_list,
  time_to_event_generator_list = ex4_event_list,
  end_fu = 10
)
ex4_data
#> # A tibble: 200 × 7
#>        id     g     x  time     f     l     y
#>     <int> <int> <int> <int> <dbl> <int> <int>
#>   1     1     0     1     1 0.173     1     0
#>   2     1     0     1     2 0.177     1     0
#>   3     1     0     1     3 0.171     1     0
#>   4     1     0     1     4 0.161     1     0
#>   5     1     0     1     5 0.163     1     0
#>   6     1     0     1     6 0.151     1     0
#>   7     1     0     1     7 0.160     1     0
#>   8     1     0     1     8 0.158     1     0
#>   9     1     0     1     9 0.179     1     0
#>  10     1     0     1    10 0.180     1     0
#>  11     2     1     0     1 0.176     0     1
#>  12     2     1     0     2 0.176     0    NA
#>  13     2     1     0     3 0.183     0    NA
#>  14     2     1     0     4 0.171     0    NA
#>  15     2     1     0     5 0.175     0    NA
ex4_data |>
  filter(time == 1) |>
  pull(l) |>
  mean()
#> [1] 0.35
ex4_data |>
  filter(time == 10) |>
  pull(l) |>
  mean()
#> [1] 0.7


# ** Example 5, realistic data set

set.seed(10956223)

ex5_fixed_list <- list(
  age = \(n) gen_ages(n = n, mean = 37, sd = 13, min = 67),
  sex = \(n) sample(c("F", "M"), size = n, replace = TRUE, prob = c(0.6, 0.4)),
  urban = \(n) rbinom(n = n, size = 1, prob = 0.8),
  region = \(n) sample(c("S", "W", "NE", "MW"), size = n, replace = TRUE, prob = c(0.38, 0.24, 0.21, 0.17))
)
ex5_varying_list <- list(
  htn = \(n, end_fu) gen_step_series(n = n, end_fu = end_fu, S1 = 0.4, t1 = 1, S2 = 0.2, t2 = 60),
  obese = \(n, end_fu) gen_step_series(n = n, end_fu = end_fu, S1 = 0.5, t1 = 1, S2 = 0.4, t2 = 60),
  copd = \(n, end_fu) gen_step_series(n = n, end_fu = end_fu, S1 = 0.8, t1 = 1, S2 = 0.75, t2 = 60),
  ulcer = \(n, end_fu) gen_step_series(n = n, end_fu = end_fu, S1 = 0.95, t1 = 1, S2 = 0.93, t2 = 60),
  unctrl_dm = \(n, end_fu) gen_time_varying_series(n = n, end_fu = end_fu, mean = 0.4)
)
ex5_event_list <- list(
  death = \(n, end_fu) gen_exp_series(n = n, end_fu = end_fu, S1 = 0.7, t1 = 60),
  disenrol = \(n, end_fu) gen_unif_series(n = n, end_fu = end_fu, S1 = 0.6, t1 = 60),
  hlth_stus = \(n, end_fu) gen_exp_series(n = n, end_fu = end_fu, S1 = 0.80, t1 = 60),
  nonadher_prota = \(n, end_fu) gen_weib2_series(n = n, end_fu = end_fu, S1 = 0.95, t1 = 12, S2 = 0.90, t2 = 60),
  nonadher_protb = \(n, end_fu) gen_weib2_series(n = n, end_fu = end_fu, S1 = 0.75, t1 = 12, S2 = 0.50, t2 = 60),
  nonadher_protc = \(n, end_fu) gen_weib2_series(n = n, end_fu = end_fu, S1 = 0.25, t1 = 12, S2 = 0.05, t2 = 60),
  ami = \(n, end_fu) gen_weib2_series(n = n, end_fu = end_fu, S1 = 0.90, t1 = 12, S2 = 0.83, t2 = 60)
)

ex5_data <- simulate_data(
  n = 1000,
  time_fixed_generator_list = ex5_fixed_list,
  time_varying_generator_list = ex5_varying_list,
  time_to_event_generator_list = ex5_event_list,
  end_fu = 60
)
ex5_data
#> # A tibble: 60,000 × 18
#>       id   age sex   urban region  time   htn obese  copd ulcer unctrl_dm
#>    <int> <dbl> <chr> <int> <chr>  <int> <int> <int> <int> <int>     <dbl>
#>  1     1  68.5 M         1 NE         1     1     0     0     0         0
#>  2     1  68.5 M         1 NE         2     1     0     0     0         0
#>  3     1  68.5 M         1 NE         3     1     0     0     0         0
#>  4     1  68.5 M         1 NE         4     1     0     0     0         0
#>  5     1  68.5 M         1 NE         5     1     0     0     0         0
#>  6     1  68.5 M         1 NE         6     1     0     0     0         0
#>  7     1  68.5 M         1 NE         7     1     0     0     0         0
#>  8     1  68.5 M         1 NE         8     1     0     0     0         0
#>  9     1  68.5 M         1 NE         9     1     0     0     0         0
#> 10     1  68.5 M         1 NE        10     1     0     0     0         0
#> # ℹ 59,990 more rows
#> # ℹ 7 more variables: death <int>, disenrol <dbl>, hlth_stus <int>,
#> #   nonadher_prota <int>, nonadher_protb <int>, nonadher_protc <int>,
#> #   ami <int>
#> # ℹ Use `print(n = ...)` to see more rows
