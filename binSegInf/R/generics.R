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

#' Generic function to form contrast vectors
#'
#' @param obj object
#' @param ... additional parameters
#'
#' @return vector of numerics
#' @export
contrast_vector <- function(obj, ...) {UseMethod("contrast_vector")}

jump_cusum <- function(obj, ...) {UseMethod("jump_cusum")}

.list_comparison <- function(obj) UseMethod(".list_comparison")

