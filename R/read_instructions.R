
#' @importFrom yaml read_yaml
#' @export read_instructions

read_instructions <- function(file) {
  
  stopifnot("Instruction file does not exist!" = base::file.exists(file))
  
  instr <- yaml::read_yaml(file = file)
  if (!"output_directory" %in% base::names(instr)) {
    warning("No output directory provided. Unable to save results.")
  }
  
  if (base::endsWith(instr$output_directory, "/")) {
    warning("Please remove the trailing '/' in the output_directory path!")
    instr$output_directory <- base::gsub("/$", "", instr$output_directory)
  }
  
  return(instr)
}


