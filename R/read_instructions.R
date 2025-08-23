#' @export read_instructions

read_instructions <- function(instr) {
  
  TF_codes <- c("TRUE" = T, "T" = T, "True" = T, "yes" = T, 
                "FALSE" = F, "F" = F, "False" = F, "no" = F)
  

  stopifnot("Please provide a scores file using the 'scores_file' argument!" = 
              "scores_file" %in% names(instr), 
            "Given scores file does not exist!" = file.exists(instr$scores_file))
  
  if (!"output_directory" %in% names(instr)) {
    warning("No output directory provided. Unable to save results.")
  }
  
  if (!"FDR" %in% names(instr) || 
      !instr$FDR %in% c("BH", "bonferroni")) {
    instr$FDR <- "BH"
  }
  
#  if (!"output_prefix" %in% names(instr) || 
#      instr$output_prefix == "") {
#    instr$output_prefix <- gsub(".csv$|.rds$", "", 
#                                basename(instr$scores_file))
#  }
  
  if (!"overwrite_output" %in% names(instr) || !instr$overwrite_output %in% TF_codes) {
    instr$overwrite_output <- T
  } else {
    instr$overwrite_output <- TF_codes[instr$overwrite_output]
  }
  

  if (!"verbose" %in% names(instr) || !instr$verbose %in% TF_codes) {
    instr$verbose <- F
  } else {
    instr$verbose <- TF_codes[as.character(instr$verbose)]
  }
  
  if (endsWith(instr$output_directory, "/")) {
    warning("Please remove the trailing '/' in the output_directory path!")
    instr$output_directory <- gsub("/$", "", instr$output_directory)
  }
  
  return(instr)
}


