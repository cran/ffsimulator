#' Connects ff_scoringhistory to past ADP rankings
#'
#' The backbone of the ffsimulator resampling process is coming up with a population of weekly outcomes for every preseason positional rank. This function creates that dataframe by connecting historical FantasyPros.com rankings to nflfastR-based scoring data, as created by `ffscrapr::ff_scoringhistory()`.
#'
#'
#' @param scoring_history a scoring history table as created by `ffscrapr::ff_scoringhistory()`
#' @param gp_model either "simple" or "none" - simple uses the average games played per season for each position/adp combination, none assumes every game is played.
#' @param pos_filter a character vector: filter the positions returned to these specific positions, default: c("QB","RB","WR","TE)
#'
#' @return a dataframe with position, rank, probability of games played, and a corresponding nested list per row of all week score outcomes.
#'
#' @examples
#' \donttest{
#' # cached data
#' scoring_history <- .ffs_cache("mfl_scoring_history.rds")
#'
#' ffs_adp_outcomes(scoring_history, gp_model = "simple")
#' ffs_adp_outcomes(scoring_history, gp_model = "none")
#' }
#'
#' @seealso `fp_rankings_history` for the included historical rankings
#' @seealso `fp_injury_table` for the historical injury table
#' @seealso `vignette("custom")` for usage details.
#'
#' @export
ffs_adp_outcomes <- function(scoring_history,
                             gp_model = "simple",
                             pos_filter = c("QB", "RB", "WR", "TE")) {
  # ASSERTIONS #
  checkmate::assert_choice(gp_model, choices = c("simple", "none"))
  checkmate::assert_character(pos_filter)
  checkmate::assert_data_frame(scoring_history)
  assert_columns(scoring_history, c("gsis_id", "team", "season", "points"))

  gsis_id <- NULL
  fantasypros_id <- NULL
  pos <- NULL
  rank <- NULL
  points <- NULL
  week <- NULL
  week_outcomes <- NULL
  player_name <- NULL

  sh <- data.table::as.data.table(scoring_history)[!is.na(gsis_id) & week <= 17,c("gsis_id", "team", "season", "points")]
  fp_rh <- data.table::as.data.table(ffsimulator::fp_rankings_history)[,-"page_pos"]
  dp_id <- data.table::as.data.table(ffscrapr::dp_playerids())[!is.na(gsis_id) & !is.na(fantasypros_id),c("fantasypros_id","gsis_id")]

  ao <- fp_rh[dp_id,on = "fantasypros_id", nomatch = 0
  ][!is.na(gsis_id) & pos %in% pos_filter
  ][sh, on = c("season","gsis_id"), nomatch = 0
  ][,list(week_outcomes = list(points), games_played = .N),
    by = c("season","pos","rank","fantasypros_id","player_name")
  ][,rank := lapply(rank, .ff_triplicate)]

  ao <- tidytable::unnest.(ao,"rank",.drop = FALSE) %>%
    .ff_apply_gp_model(gp_model)

  ao <- ao[
    ,list(week_outcomes = list(c(unlist(week_outcomes))),
          player_name = list(player_name),
          fantasypros_id = list(fantasypros_id)
          ),
    by = c("pos","rank","prob_gp")
    ][order(pos,rank)][!is.na(fantasypros_id)]

  return(ao)
}

#' Applies various injury models to adp outcomes
#'
#' @keywords internal
#' @return same adp outcomes dataframe but with a prob_gp column
.ff_apply_gp_model <- function(adp_outcomes, model_type) {
  if (model_type == "none") {
    adp_outcomes$prob_gp <- 1
  }

  if (model_type == "simple") {
    adp_outcomes <- adp_outcomes[data.table::as.data.table(ffsimulator::fp_injury_table),on = c("pos","rank")]
  }

  adp_outcomes
}

.ff_triplicate <- function(.x){
  c(ifelse(.x-1==0,.x,.x-1),.x,.x+1)
}
