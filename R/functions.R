#' @export

read_instructions <- function(yaml_fpath) {
  
  stopifnot("Instruction file does not exist!" = file.exists(yaml_fpath))
  
  instr <- read_yaml(file = yaml_fpath)
  
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




#usable_for_limma <- \(.x) {
#  .usable <- is.array(.x) & sum(is.na(.x)) / length(.x) <= 0.9
#  return(.usable)
#}























collect_all_layer_configurations <- function(GI_data, 
                                             .to_use = c("tech_rep", "bio_rep", "guide_pair"), 
                                             make_pos_agnostic = F) {
  
  .collapsable_layers <- intersect(c("tech_rep", "bio_rep", "guide_pair"), colnames(GI_data))
  .to_use = intersect(.to_use, colnames(GI_data))
  
  output <- list()
  
  for (.use in .to_use) {
    
    output[[str_c("default_", .use, "_used")]] <- GIScores(
      GI_data, 
      block_layer = .use, 
      pos_agnostic = make_pos_agnostic)
  }
  
  for (.i in seq_len(length(.collapsable_layers)-1)) {
    for (.to_collapse in combn(.collapsable_layers, .i, simplify = F)) {
      for (.use in setdiff(.collapsable_layers, .to_collapse)) {
        .n <- str_c(str_c(.to_collapse, collapse = "_"), "_collapsed_", .use, "_used")
        print(.n)
        output[[.n]] <- GIScores(GI_data, 
                                 collapse_layers = .to_collapse, 
                                 block_layer = .use, 
                                 pos_agnostic = make_pos_agnostic)
      }}}
  return(output)
}




find_optimal_configuration <- function(GI_list, verbose = T) {
  
  print(str(map(GI_list, ~ .x@dupCorrelation)))
  
  .x <- GI_list |> map(~ .x@dupCorrelation) |> map_dbl(mean, na.rm = T)
  
  .d <- data.table(config = names(GI_list), 
                   dcor = GI_list |> map(~ .x@dupCorrelation) |> map_dbl(mean, na.rm = T) |> round(3))
  
  
  if (any(.x >= 0)) {.x <- keep(.x, .x >= 0)} # if greater or = 0 found, remove <= 0
  if (any(.x > 0)) {.x <- .x |> keep(.x > 0)} # if any greater 0 found, remove = 0
  if (length(.x) >= 1) {.x <- .x[which.min(.x)]} # if more than one option remains, choose weakest correlation
  
  .d[, kept := fcase(config %in% names(.x), "*", default = "")]
  
  for (.n in names(GI_list)) {
    GI_list[[.n]]@metadata$dupcor_data <- .d
    
    GI_list[[.n]]@metadata$dupcor_plot <- ggplot(
      data = .d, mapping = aes(dcor, config)) + 
      theme_light() + 
      geom_col(aes(fill = kept)) + 
      scale_color_manual(values = c("TRUE" = "seagreen", "FALSE" = "NA")) + 
      geom_vline(xintercept = c(0, 0.25), linetype = "dashed", linewidth = 1) +
      labs(#caption = "It is recommended to choose a configuration with most values between 0 and 0.25.", 
        x = "Duplicate correlation", 
        y = "Limma configuration")
  }
  #if (verbose) {
  print(.d)
  #  }
  
  return(GI_list[names(.x)])
}





