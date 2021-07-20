## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
options(dplyr.summarise.inform = FALSE)

## ----setup, message = FALSE---------------------------------------------------
library(ffsimulator)
library(ggplot2)

## ----eval = FALSE-------------------------------------------------------------
#  foureight_conn <- mfl_connect(2021, 22627)
#  
#  foureight_sim <- ff_simulate(conn = foureight_conn, n_seasons = 25, n_weeks = 14)
#  foureight_sim

## ----echo = FALSE, message = FALSE--------------------------------------------
foureight_sim <- readRDS(system.file("cache/foureight_sim.rds", package = "ffsimulator"))
foureight_sim

## -----------------------------------------------------------------------------
autoplot(foureight_sim) # defaults to type = "wins"
autoplot(foureight_sim, type = "rank")
autoplot(foureight_sim, type = "points")

## -----------------------------------------------------------------------------
foureight_sim$summary_simulation

## -----------------------------------------------------------------------------
foureight_sim$summary_season

## -----------------------------------------------------------------------------
foureight_sim$summary_week

## -----------------------------------------------------------------------------
foureight_sim$roster_scores

## -----------------------------------------------------------------------------
foureight_sim$league_info

