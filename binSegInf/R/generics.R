#' Checks S3 object validity
#' 
#' Generic function that fails noisily with a stop message if the object is 
#' invalid. Otherwise, nothing happens.
#' 
#' @param  obj  The object to check
#' @return void
#' @export
isValid <- function(obj) UseMethod("isValid")

form_polyhedra <- function(obj, ...) UseMethod("form_polyhedra")

.list_comparison <- function(obj) UseMethod(".list_comparison")

