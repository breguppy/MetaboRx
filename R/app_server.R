#' @keywords internal
#' @noRd

app_server <- function(input, output, session) {
  session$onSessionEnded(function() { stopApp() })
  
  import <- mod_import_server("import")
  
  base_data <- reactive({
    .merge_lists(list(
      cleaned  = .get_or_null(import$cleaned),
      filtered = .get_or_null(import$filtered)
    ))
  })
  
  correct <- mod_correct_server("correct",
                                data   = base_data,
                                params = import$params)
  
  base_params <- reactive( .merge_lists(import$params(), correct$params()) )
  combined_data <- reactive({
     .merge_lists(
       base_data(),
       list(
         imputed            = .get_or_null(correct$imputed),
         corrected          = .get_or_null(correct$corrected),
         filtered_corrected = .get_or_null(correct$filtered_corrected),
         transformed        = .get_or_null(correct$transformed),
         all_corr            = .get_or_null(correct$all_corr)
       )
     )
   })
   
   viz <- mod_visualize_server("viz",
                               data   = combined_data,
                               params = base_params)
   
  combined_params <- reactive( .merge_lists(base_params(), viz$params()) )
 
  mod_export_server("export",
                   data   = combined_data,
                   params = combined_params)
}
