## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  eval = FALSE
)
options(dplyr.summarise.inform = FALSE)

## ----setup, message = FALSE---------------------------------------------------
#  library(ffsimulator)
#  library(ffscrapr)
#  library(dplyr)
#  library(ggplot2)

## ----eval = FALSE-------------------------------------------------------------
#  scoring_history <- ffscrapr::ff_scoringhistory(conn, seasons = 2012:2020)

## ----eval = FALSE-------------------------------------------------------------
#  latest_rankings <- ffs_latest_rankings()

## ----eval = FALSE-------------------------------------------------------------
#  rosters <- ffs_rosters(conn)

## ----eval = FALSE-------------------------------------------------------------
#  lineup_constraints <- ffscrapr::ff_starter_positions(conn)

## ----eval = FALSE-------------------------------------------------------------
#  league_info <- ffscrapr::ff_league(conn)

## ----eval = FALSE-------------------------------------------------------------
#  adp_outcomes <- ffs_adp_outcomes(
#      scoring_history = scoring_history,
#      injury_model = "simple" # or "none"
#    )

## ----eval = FALSE-------------------------------------------------------------
#  projected_scores <- ffs_generate_projections(
#      adp_outcomes = adp_outcomes,
#      latest_rankings = latest_rankings,
#      n_seasons = 100, # number of seasons
#      n_weeks = 14, # weeks per season
#      rosters = rosters # optional, reduces the sample to just rostered players
#    )

## ----eval = FALSE-------------------------------------------------------------
#  roster_scores <- ffs_score_rosters(
#      projected_scores = projected_scores,
#      rosters = rosters
#    )

## ----eval = FALSE-------------------------------------------------------------
#  optimal_scores <- ffs_optimise_lineups(
#      roster_scores = roster_scores,
#      lineup_constraints = lineup_constraints,
#      lineup_efficiency_mean = 0.775,
#      lineup_efficiency_sd = 0.05,
#      best_ball = FALSE, # or TRUE
#      parallel = FALSE # or TRUE
#    )

## ----eval = FALSE-------------------------------------------------------------
#    schedules <- ffs_build_schedules(
#      n_teams = length(unique(rosters$franchise_id)),
#      n_seasons = n_seasons,
#      n_weeks = n_weeks,
#      seed = NULL
#    )

## ----eval = FALSE-------------------------------------------------------------
#  summary_week <- ffs_summarise_week(optimal_scores, schedules)
#  summary_season <- ffs_summarise_season(summary_week)
#  summary_simulation <- ffs_summarise_simulation(summary_season)

## ----eval = FALSE-------------------------------------------------------------
#  options(ffscrapr.cache = "filesystem")
#  library(ffsimulator)
#  library(ffscrapr)
#  library(tidyverse)
#  library(furrr)
#  library(tictoc)
#  
#  plan(multisession)
#  set.seed(613)

## ----eval = FALSE-------------------------------------------------------------
#  conn <- mfl_connect(2021, 47747) # a random SFB league to grab league info from
#  
#  league_info <- ffscrapr::ff_league(conn)
#  
#  scoring_history <- ffscrapr::ff_scoringhistory(conn, 2012:2020)
#  
#  adp_outcomes <- ffs_adp_outcomes(scoring_history = scoring_history, injury_model = "simple")
#  
#  latest_rankings <- ffs_latest_rankings()
#  lineup_constraints <- ffscrapr::ff_starter_positions(conn)

## ----eval = FALSE-------------------------------------------------------------
#  conn2 <- mfl_connect(2021)
#  
#  leagues <- mfl_getendpoint(conn2, "leagueSearch", SEARCH = "#SFB11") %>%
#    pluck("content","leagues","league") %>%
#    tibble() %>%
#    unnest_wider(1) %>%
#    filter(str_detect(name,"Mock|Copy|Satellite|Template",negate = TRUE))
#  
#  
#  get_rosters <- function(league_id){
#    mfl_connect(2021, league_id) %>%
#      ffs_rosters()
#  }
#  get_franchises <- function(league_id){
#    mfl_connect(2021, league_id) %>%
#      ff_franchises()
#  }
#  
#  rosters_raw <- leagues %>%
#    select(-homeURL) %>%
#    mutate(
#      rosters = future_map(id, get_rosters),
#      franchises = future_map(id, get_franchises)
#    )
#  
#  franchises <- rosters_raw %>%
#    select(league_id = id, franchises) %>%
#    unnest(franchises) %>%
#    select(league_id, franchise_id, division_name)
#  
#  rosters <- rosters_raw %>%
#    select(rosters) %>%
#    unnest(rosters) %>%
#    left_join(franchises,by = c("league_id","franchise_id"))

## ----eval = FALSE-------------------------------------------------------------
#  n_seasons <- 100
#  n_weeks <- 13
#  projected_scores <- ffs_generate_projections(adp_outcomes = adp_outcomes,
#                                               latest_rankings = latest_rankings,
#                                               n_seasons = n_seasons,
#                                               n_weeks = n_weeks,
#                                               rosters = rosters)
#  
#  tictoc::tic(glue::glue("ffs_score_rosters {Sys.time()}"))
#  roster_scores <- ffs_score_rosters(projected_scores, rosters)
#  tictoc::toc()
#  
#  tictoc::tic("ffs_optimize_lineups {Sys.time()}")
#  optimal_scores <- ffs_optimize_lineups(
#    roster_scores = roster_scores,
#    lineup_constraints = lineup_constraints,
#    best_ball = FALSE,
#    parallel = TRUE)
#  tictoc::toc()

## ----eval = FALSE-------------------------------------------------------------
#  schedules <- ffs_build_schedules(n_teams =
#                                     rosters %>%
#                                      dplyr::distinct(
#                                        .data$league_id,
#                                        .data$franchise_id) %>%
#                                      nrow(),
#                                   n_seasons = n_seasons,
#                                   n_weeks = n_weeks)
#  
#  summary_week <- ffs_summarise_week(optimal_scores, schedules)
#  summary_season <- ffs_summarise_season(summary_week)
#  summary_simulation <- ffs_summarise_simulation(summary_season)

