#' Checks S3 object validity
#' 
#' Generic function that fails noisily with a stop message if the object is 
#' invalid. Otherwise, nothing happens.
#' 
#' @param  obj  The object to check
#' @return void
#' @export
isValid <- function(obj) UseMethod("isValid")

#' Generic function to generate polyhedra
#'
#' @param obj object
#' @param ... additional parameters
#'
#' @return An object of class polyhedra
#' @export
form_polyhedra <- function(obj, ...) {UseMethod("form_polyhedra")}

#' Generic function to count the number of jumps
#'
#' @param obj object
#' @param ... additional parameters
#'
#' @return vector of jumps
#' @export
get_jumps <- function(obj, ...) {UseMethod("get_jumps")}

#' Generic function to form contrast vectors
#'
#' @param obj object
#' @param ... additional parameters
#'
#' @return vector of numerics
#' @export
contrast_vector <- function(obj, ...) {UseMethod("contrast_vector")}

get_jump_cusum <- function(obj, ...) {UseMethod("get_jump_cusum")}

.list_comparison <- function(obj) UseMethod(".list_comparison")

