#' @export full_run


full_run <- function(yaml_fpath, return_output = F) {
  
  stopifnot("Instruction file does not exist!" = file.exists(yaml_fpath))
  
  instr <- read_instructions(yaml::read_yaml(file = yaml_fpath))

  #for (.pkg in c("BiocManager")) {if (!require(.pkg, quietly = T)) {utils::install.packages(.pkg)}}
  #if (!require("limma", quietly = T)) {BiocManager::install("limma")}
  #if (F) {devtools::load_all("D:/Promotion_projects_github/CeRberus/")}
  
  print(str(instr))
  
  .data <- switch(tools::file_ext(instr$scores_file),
                  "csv" = data.table::fread(instr$scores_file), 
                  "rds" = readRDS(instr$scores_file))
  
  .data <- GIScores(.data)
  
  .data <- new_limma_object(.data)
  
  if (instr$verbose) {
    print("(1/5) Data import successful.")
    print("Deciding on block structure. This may take a while.")
  }
  
  
  
  while(is.null(attr(.data, "block_layer")) & 
        ncol(attr(.data, "replicate_layers")) > 1) {
    
    # This should find the block with the highest in-block correlation. 
    # If no positive block is found, the lowest correlating layer gets collapsed, and the process is repeated. 
    
    .dc <- purrr::map_dbl(compute_dupcor_values(.data, sample_query = 20), median, na.rm = T)
    
    .highest <- .dc[which.max(.dc)]
    
    if (.highest >= 0) {
      data.table::setattr(.data, "block_layer", names(.highest))
    } else {
      .data <- collapse_layer(.data, names(which.min(.dc)))
    }
  }
  
  
  if (instr$verbose) print("(2/5) Optimal block structure identified. Estimating duplicate correlation.")
  
  .data <- list(data = .data, 
                dupcor = compute_dupcor_values(.data))
  
  if (instr$verbose) print("(3/5) Computed duplicate correlation values.")
  
  .data$GI_scores <- compute_GIs(GI_object = .data$data, 
                                 dupcor = .data$dupcor, 
                                 FDR_method = instr$FDR)
  
  
  #print(str(.data$GI_scores))
  
  if (instr$verbose) print("(4/5) Computed GI scores.")
  
  
  
  print(str(.data$GI_object))
  
  
  if ("output_directory" %in% names(instr) & instr$overwrite_output) {
    
    export_GIs(
      GI_object = .data$GI_scores, 
      dupcor_object = .data$dupcor, 
      directory = file.path(instr$output_directory)
    )
    
    .log <- c(list(instructions = instr), 
              list(file1 = .data$data, 
                   file3 = .data$dupcor) |> map(~ attributes(.x)[!names(attributes(.x)) %in% c("dim", "dimnames"#, "dim_description", "replicate_layers"
                   )]), 
              screen_design = list(
                queried = attr(.data$data, "queried"), 
                collapsed_layers = if (is.null(attr(.data$GI_scores, "collapsed_layers"))) {"none"} else {attr(.data$GI_scores, "collapsed_layers")}, 
                block_layer = if (is.null(attr(.data$GI_scores, "block_layer"))) {"none"} else {attr(.data$GI_scores, "block_layer")}
              ))
    
    #print(str(.log))
    
    yaml::write_yaml(.log, file.path(instr$output_directory, "limma.log"))
    
    
    if (instr$verbose) print("(5/5) Exported results. Limma run successful!")
  } else {
    if (instr$verbose) print("(5/5) Limma run successful. Output not saved.")
  }
  
  if (return_output) {
    return(list(data = .data, log = .log))
  } else {
    return(NULL)
    }
}