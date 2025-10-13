#' @export full_run


full_run <- function(yaml_fpath, return_output = F) {
  
  stopifnot("Instruction file does not exist!" = file.exists(yaml_fpath))
  
  instr <- read_instructions(yaml::read_yaml(file = yaml_fpath))
  
  stopifnot("Output directory required!" = !is.null(instr$output_directory))
  
  dir.create(instr$output_directory, showWarnings = F, recursive = T)
  
  #for (.pkg in c("BiocManager")) {if (!require(.pkg, quietly = T)) {utils::install.packages(.pkg)}}
  #if (!require("limma", quietly = T)) {BiocManager::install("limma")}
  #if (F) {devtools::load_all("D:/Promotion_projects_github/CeRberus/")}
  
  print(str(instr))
  
  .data <- switch(tools::file_ext(instr$scores_file),
                  "csv" = data.table::fread(instr$scores_file), 
                  "rds" = readRDS(instr$scores_file))
  
  .data <- GIScores(.data)
  
  if (instr$verbose) {
    print("(1/5) Data import successful.")
    print("Deciding on block structure. This may take a while.")
  }
  
  .data <- block_decision_heuristics(.data)
  
  if (instr$verbose) print("(2/5) Optimal block structure identified. Estimating duplicate correlation.")
  
  .data <- compute_dupcor_values(.data)
  
  if (instr$verbose) print("(3/5) Computed duplicate correlation values.")
  
  .data <- compute_GIs(.data, FDR_method = instr$FDR)
  
  
  if (instr$verbose) print("(4/5) Computed GI scores.")
  
  
  if ("output_directory" %in% names(instr) & instr$overwrite_output) {
    
    .output <- list()
    
    .output[[str_glue("GI_scores")]] <- GI_df(.data, .block = blocks(.data)$chosen)
    .output[[str_glue("duplicate_correlation")]] <- dupCorrelation_df(.data, .block = blocks(.data)$chosen)
    
    
    
    .output |>
      iwalk(~ {data.table::fwrite(
        x = .x, 
        file = file.path(instr$output_directory, paste0(.y, ".csv")))
      })
    
    .log <- create_log(.data)
    writeLines(.log, file.path(instr$output_directory, "limma.log"))
    
    if (instr$verbose) print("(5/5) Exported results. Limma run successful!")
  } else {
    if (instr$verbose) print("(5/5) Limma run successful. Output not saved.")
  }
  
  if (return_output) {return(.data)} else {return(NULL)}
}