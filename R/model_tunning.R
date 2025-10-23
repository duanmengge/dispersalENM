# feature combination
# regularization multiplier
#' model_tuning function
#'
#' @title Select the optimal parameter among the alternative parameters of
#' MaxEnt model
#' @description In the Maxent model, regularization multiplier and feature
#' combination are two important parameters. By model tuning, the optimal
#' combination is calculated after a given set of alternative parameter
#' combination.
#' @usage
#' model_tuning(occs,
#'              env,
#'              nbg = 10000,
#'              f_c = c("L", "Q", "P"),
#'              r_m = c(1, 1, 3),
#'              train_ratio = 0.7)
#' @param occs data.frame. species occurrence data.
#' @param env SpatRaster. environmental data.
#' @param nbg numeric. The quantity of background data. Default is 10000.
#' @param f_c feature combination. The option including "L"(linear),
#' "Q"(quadratic), "P"(product), "T"(threshold), "H"(hinge). Combine these
#' selected options and optimize the model. the default is c("L","Q","P")
#' @param r_m regularization multiplier. The format like default c(1, 1, 3).
#' The first number is the starting value, the third number is the termination
#' value, and the second number is step.
#' @param train_ratio 0-1. Proportion of training data
#'
#' @return a data frame
#' @export
#'
#' @examples
#' occourence <- readRDS(system.file("extdata", "occourence.rda",
#'                                  package = "dispersalENM"))
#' trainenv <-  terra::rast(list.files(
#' system.file("extdata", "env/trainenv", package = "dispersalENM"),
#' pattern = "tif", full.names = TRUE))
#' #model_tuning
#' model_tuning_results <- model_tuning(occourence, trainenv, nbg=10000,
#'                                        f_c=c("L","Q"), r_m=c(1, 1, 3))
model_tuning <- function(occs,
                          env,
                          nbg = 10000,
                          f_c = c("L", "Q", "P"),
                          r_m = c(1, 1, 3),
                          train_ratio = 0.7,
                          partitions = "block") {
  if (partitions == "testing") {
    train <- sample(nrow(occs), train_ratio * nrow(occs))
    occ_train <- occs[train, ]
    occ_test <- occs[-train, ]
  } else {
    occ_train <- occs
    occ_test <- NULL
  }
  get_combinations <- function(feature) {
    if (length(feature) == 1) {
      return(as.character(feature))
    } else {
      prev_combinations <- get_combinations(feature[-length(feature)])
      new_combinations <- c(
        prev_combinations,
        paste(prev_combinations, feature[length(feature)], sep = "")
      )
      return(unique(c(as.character(feature), new_combinations)))
    }
  }
  fc_combinations <- get_combinations(f_c)
  print(fc_combinations)
  rm_combinations <- seq(r_m[1], r_m[3], by = r_m[2])
  print(rm_combinations)
  tune_args <- list(fc = fc_combinations, rm = rm_combinations)
  other_settings <- list(
    abs.auc.diff = FALSE, pred.type = "logistic",
    validation.bg = "partition"
  )
  partition_settings <- list(orientation = "lon_lat")
  e <- ENMeval::ENMevaluate(
    occs = occ_train, occs.testing = occ_test, env, n.bg = nbg,
    tune.args = tune_args, partitions = partitions,
    other.settings = other_settings, partition.settings = partition_settings,
    algorithm = "maxnet", overlap = FALSE, doClamp = FALSE
  )
  results <- e@results
  return(results)
}
