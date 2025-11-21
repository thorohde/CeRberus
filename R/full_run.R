#' @export


full_run <- function(yaml_fpath, return_output = F) {
  
  stopifnot("Instruction file does not exist!" = file.exists(yaml_fpath))
  
  instr <- read_instructions(yaml::read_yaml(file = yaml_fpath))
  
  stopifnot("Output directory required!" = !is.null(instr$output_directory))
  
  dir.create(instr$output_directory, showWarnings = F, recursive = T)
  
  #for (.pkg in c("BiocManager")) {if (!require(.pkg, quietly = T)) {utils::install.packages(.pkg)}}
  #if (!require("limma", quietly = T)) {BiocManager::install("limma")}
  #if (F) {devtools::load_all("D:/Promotion_projects_github/CeRberus/")}
  
  if (instr$verbose) {print(str(instr, give.attr = F))}
  
  .data <- switch(tools::file_ext(instr$scores_file),
                  "csv" = data.table::fread(instr$scores_file), 
                  "rds" = readRDS(instr$scores_file))
  
  
  if (instr$verbose) {
    print("(1/5) Data import successful.")
    
  }
  
  if (instr$verbose) {print("Deciding on block structure. This may take a while.")}  
  
  .data <- collect_all_layer_configurations(.data)
  
  
  
  library(ggplot2)
  
  .plt_data <- data.frame(config = names(.data), dupcor = map_dbl(.data, dupCorrelation))
  
  print(.plt_data)
  
  if (!"keep_all_configurations" %in% names(instr) || (!instr$keep_all_configurations)) {
    .data <- find_optimal_configuration(.data, verbose = instr$verbose)
  }
  
  .plt_data$kept <- .plt_data$config %in% names(.data)
  
  
  dupcor_plot <- \(plt_data, fname = NULL) {
    
    plt <- ggplot(data = plt_data, 
                  mapping = aes(dupcor, config)) + 
      theme_light() + 
      geom_col(aes(fill = kept)) + 
      scale_color_manual(values = c("TRUE" = "seagreen", "FALSE" = "NA")) + 
      geom_vline(xintercept = c(0, 0.25), linetype = "dashed", linewidth = 1) +
      labs(#caption = "It is recommended to choose a configuration with most values between 0 and 0.25.", 
        x = "Duplicate correlation", 
        y = "Limma configuration")
    
    if (!is.null(fname)) {
      ggplot2::ggsave(filename = fname, 
                      plot = plt, width = 8, 
                      height = 5, 
                      dpi = 300)}
  }
  
  
  dupcor_plot(.plt_data, fname = file.path(instr$output_directory, "duplicateCorrelationPlot.png"))
  
  
  if (instr$verbose) print("(2/5) Optimal block structure identified. Estimating duplicate correlation.")
  
  .data <- map(.data, compute_GIs, FDR_method = instr$FDR)
  
  
  if (any(map_chr(.data, screenType) == "multiplex.symmetric.position.agnostic")) {
    if (instr$verbose) {message("Computing symmetric GIs. This might take a while.")}
    
    .data <- modify_if(.data, grepl("position.agnostic", map(.data, screenType)), compute_symmetric_GIs)} # FDR mmethod
  
  
  if (instr$verbose) print("(4/5) Computed GI scores.")
  
  if ("output_directory" %in% names(instr) & instr$overwrite_output) {
    
    .output <- list()
    
    print(names(.data))
    
    .output$duplicate_correlation <- .plt_data
    
    for (.n in names(.data)) {
      .output[[paste0("GI_scores_", .n)]] <- GI_df(.data[[.n]])

      if (screenType(.data[[.n]]) == "multiplex.symmetric.position.agnostic") {
        .output[[paste0("GI_scores_position_agnostic_", .n)]] <- symmetricGI_df(.data[[.n]])
      }
    }
    
    .output |>
      iwalk(~ {data.table::fwrite(
        x = .x, 
        file = file.path(instr$output_directory, paste0(.y, ".csv")))
      })
    
    #.log <- create_log(.data)
    
    #writeLines(.log, file.path(instr$output_directory, "limma.log"))
    
    
    
    
    if (instr$verbose) print("(5/5) Exported results. Limma run successful!")
  } else {
    if (instr$verbose) print("(5/5) Limma run successful. Output not saved.")
  }
  
  if (return_output) {return(.data)} else {return(NULL)}
}