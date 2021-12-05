#' Append file name to path
#'
#' Path is in enclosing environment!
#' Helper for running scripts.
#'
#' @param path file path to append to \code{path}
#'
#' @return New file path
#' @export
#'
#' @examples
#' path = "../local_data_ukko2/"
#' pap = pap_factory(path)
#' pap("myfile.whatever")
pap_factory <- function(path){
  function(str){
    paste0(path, str)
  }
}
