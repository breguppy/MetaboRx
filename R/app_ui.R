#' app_ui.R
#' @importFrom bslib card navset_tab nav_panel layout_sidebar sidebar bs_theme
#' @importFrom shiny tags icon fluidPage
#' @keywords internal
#' @noRd

app_ui <- function() {
  fluidPage(
    theme = bslib::bs_theme(preset = "cosmo"),
    tags$head(
      tags$style(
        HTML(
          "
      /* default state of tabs */
      .nav-tabs > li > a {
        background-color: #2780E3;   /* unselected tab bg */
        color: #FFFFFF;              /* unselected tab text */
        border-radius: 1;
      }

      .nav-tabs > li > a:hover {
        background-color: #e0e0e0;
        color: #2c3e50;
      }

      /* active tab */
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        background-color: #1E88E5;   /* selected tab bg */
        color: #ffffff;              /* selected tab text */
        border-color: #1E88E5;
      }

      .popover.popover-responsive {
    max-width: min(90vw, 900px);
    width: min(90vw, 900px);
  }

  .popover.popover-responsive .popover-body {
    max-height: 80vh;
    overflow-y: auto;
  }

  .popover.popover-responsive img {
    max-width: 100%;
    height: auto;
  }
    ")
      )
    ), shinyjs::useShinyjs(), titlePanel("MetaboRx: Quality Assessment and Preprocessing for Metabolomics Data"), bslib::navset_tab(
      id = "main_steps", mod_import_ui("import"), mod_correct_ui("correct"), mod_visualize_ui("viz"), mod_export_ui("export")
    )
  )
}
