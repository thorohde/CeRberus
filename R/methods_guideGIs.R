

setMethod("collapse_replicates", 
          signature = signature(.x = "GuideGI"), 
          function(.x) {
            
          #  print("Before: ")
          #  print(str(.x))
            
            if (l(.x@collapse) > 0) {
              stopifnot(all(.x@collapse %in% .x@replicates))
              
              .margin <- setdiff(1:length(dim(.x@data)), 
                                 which(c(.x@space, .x@replicates) %in% .x@collapse))
              
              .x@data <- apply(X = .x@data, 
                               MARGIN = .margin, 
                               FUN = mean, na.rm = T)
              
              .x@replicates <- setdiff(.x@replicates, .x@collapse)
            }
            
            
         #   print("After: ")
         #   print(str(.x))
            
            return(.x)
          })




setMethod("flatten_guideGIs", 
          signature = signature(.x = "GuideGI"), 
          function(.x) {
            
            .x@use_blocks <- .x@block_layer != "" & 
                length(.x@block_layer) > 0 & 
                length(setdiff(.x@replicates, .x@block_layer)) != 0

            .f <- as.formula(paste0(paste0(.x@space, collapse = " ~ "), " ~ ", 
                                    paste0(.x@replicates, collapse = " + ")))
            
            .x@data <- .x@data |> 
              flatten_array(c(.x@space, .x@replicates)) |>
              acast(.f)
            
            .x@block_description <- dimnames(.x@data)[[l(.x@space)+1]]
            
            if (length(.x@blocks) != 1 && !identical(.x@blocks, "none")) {
              
              .x@blocks <- c(guide_pair = "(g\\d+)", 
                             bio_rep = "(b\\d+)", 
                             tech_rep = "(t\\d+)")[[.x@block_layer]]
              
              .x@blocks <- as.character(factor(str_match(.x@block_description, .x@blocks)[,2]))
            }
            return(.x)
          })




setMethod("compute_dupCorrelation", 
          signature = signature(.x = "GuideGI"), 
          function(.x) {
            if (length(.x@space) == 1) {
              
              if (.x@use_blocks) {
                output <- suppressWarnings(
                  limma::duplicateCorrelation(
                    object = .x@data, 
                    block = .x@blocks))
                
              } else {
                output <- suppressWarnings(
                  limma::duplicateCorrelation(object = .x@data))
              }
              output <- output$consensus.correlation
            }
            
            if (length(.x@space) == 2) {
              
              output <- set_names(rownames(.x@data)) |>
                map(\(.g) .x@data[.g,,])
              
              if (.x@use_blocks) {
                output <- output |>
                  map(\(.g) suppressWarnings(
                    limma::duplicateCorrelation(
                      object = .g, 
                      block = .x@blocks)))
              } else {
                output <- output |>
                  map(\(.g) suppressWarnings(
                    limma::duplicateCorrelation(object = .g, 
                                                ndups = 1 # required 
                                                )))
                
              }
              output <- output |> map_dbl("consensus.correlation")
            }
            return(output)
          })
















