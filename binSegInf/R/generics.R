#' Checks S3 object validity
#' 
#' Generic function that fails noisily with a stop message if the object is 
#' invalid. Otherwise, nothing happens.
#' 
#' @param  obj  The object to check
#' @return void
#' @export
is_valid <- function(obj) UseMethod("is_valid")

#' Generic function to generate polyhedra
#'
#' @param obj object
#' @param ... additional parameters
#'
#' @return An object of class polyhedra
#' @export
polyhedra <- function(obj, ...) {UseMethod("polyhedra")}

#' Generic function to count the number of jumps
#'
#' @param obj object
#' @param ... additional parameters
#'
#' @return vector of jumps
#' @export
jumps <- function(obj, ...) {UseMethod("jumps")}

#' Generic function to count the number of jump cusum
#'
#' @param obj  object
#' @param ... additional parameters
#'
#' @return vector of cusums
#' @export
jump_cusum <- function(obj, ...) {UseMethod("jump_cusum")}

#' Generic function to count the number of jump lambdas
#'
#' @param obj  object
#' @param ... additional parameters
#'
#' @return vector of lambdas
#' @export
jump_lambda <- function(obj, ...) {UseMethod("jump_lambda")}

.list_comparison <- function(obj) {UseMethod(".list_comparison")}
.get_length <- function(obj) {UseMethod(".get_length")}
