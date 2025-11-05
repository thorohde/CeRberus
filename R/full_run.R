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
  
  .data <- collect_all_layer_configurations(.data)
  
  if (instr$verbose) {
    print("(1/5) Data import successful.")
    print("Deciding on block structure. This may take a while.")
  }
  
  .data <- map(.data, compute_dupcor_values, sample_query = 50, .progress = T)
  
  if (!"keep_all_configurations" %in% names(instr) || 
      (!instr$keep_all_configurations)) {
    .data <- find_optimal_configuration(.data)
  }
  
  if (instr$verbose) print("(2/5) Optimal block structure identified. Estimating duplicate correlation.")
  
  .data <- map(.data, compute_dupcor_values, .progress = T)
  
  if (instr$verbose) print("(3/5) Computed duplicate correlation values.")
  
  .data <- map(.data, compute_GIs, FDR_method = instr$FDR)
  
  if (instr$verbose) print("(4/5) Computed GI scores.")
  
  
  if ("output_directory" %in% names(instr) & instr$overwrite_output) {
    
    .output <- list()
    
    print(names(.data))
    
    for (.n in names(.data)) {
      .output[[paste0("GI_scores_", .n)]] <- GI_df(.data[[.n]])
      .output[[paste0("duplicate_correlation_", .n)]] <- dupCorrelation_df(.data[[.n]])
      
    }
    
    .output |>
      iwalk(~ {data.table::fwrite(
        x = .x, 
        file = file.path(instr$output_directory, paste0(.y, ".csv")))
      })
    
    #.log <- create_log(.data)
    
    #writeLines(.log, file.path(instr$output_directory, "limma.log"))
    
    library(ggplot2)
    
    .plt <- .data |> 
      imap(~ data.table(config = .y, dupCorrelation_df(.x))) |>
      rbindlist() |>
      ggplot(aes(dupcor, config)) + theme_light() + 
      geom_boxplot() + 
      geom_jitter(size = 0.5) + 
      geom_vline(xintercept = c(0, 0.25), linetype = "dashed", linewidth = 1) +
      labs(caption = "It is recommended to choose a configuration with most values between 0 and 0.25.", 
           x = "Duplicate correlation", 
           y = "Limma configuration")
    
    ggplot2::ggsave(filename = file.path(instr$output_directory, "duplicateCorrelationPlot.png"), 
                    plot = .plt, width = 8, 
                    height = 5, 
                    dpi = 300)
    
    
    
    
    if (instr$verbose) print("(5/5) Exported results. Limma run successful!")
  } else {
    if (instr$verbose) print("(5/5) Limma run successful. Output not saved.")
  }
  
  if (return_output) {return(.data)} else {return(NULL)}
}