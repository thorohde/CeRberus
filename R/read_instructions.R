
#' @importFrom yaml read_yaml
#' @export read_instructions

read_instructions <- function(file) {
  
  stopifnot("Instruction file does not exist!" = file.exists(file))
  
  instr <- yaml::read_yaml(file = file)
  
  stopifnot("Please provide a scores file using the 'scores_file' argument!" = 
              "scores_file" %in% names(instr), 
            "Given scores file does not exist!" = file.exists(instr$scores_file))
  
  
  if (!"output_directory" %in% names(instr)) {
    warning("No output directory provided. Unable to save results.")
  }
  
  if (!"output_prefix" %in% names(instr) || 
      instr$output_prefix == "") {
    instr$output_prefix <- gsub(".csv$|.rds$", "", 
                                basename(instr$scores_file))
  }

  if (!"overwrite_output" %in% names(instr) || 
      !instr$overwrite_output %in% c("TRUE", "T", "FALSE", "F")) {
    instr$overwrite_output <- T
  }
  
    
  if (endsWith(instr$output_directory, "/")) {
    warning("Please remove the trailing '/' in the output_directory path!")
    instr$output_directory <- gsub("/$", "", instr$output_directory)
  }
  
  return(instr)
}


