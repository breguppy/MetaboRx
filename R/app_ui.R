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
          /* ================================================================
             Normal screen styling
             ================================================================ */

          /* Default state of tabs */
          .nav-tabs > li > a {
            background-color: #2780E3;
            color: #FFFFFF;
            border-radius: 1;
          }

          .nav-tabs > li > a:hover {
            background-color: #e0e0e0;
            color: #2c3e50;
          }

          /* Active tab */
          .nav-tabs > li.active > a,
          .nav-tabs > li.active > a:focus,
          .nav-tabs > li.active > a:hover {
            background-color: #1E88E5;
            color: #ffffff;
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

          /* ================================================================
             PDF and printer styling
             Activated automatically by Ctrl+P
             ================================================================ */

          @media print {
            @page {
              size: letter landscape;
              margin: 0.5in;
            }

            /*
             * All text within a manuscript-summary container is printed
             * at exactly 12 points.
             *
             * !important overrides:
             * - the 1.4em value in metric_card()
             * - the 0.85rem value in correlation badges
             * - Bootstrap's .small class
             * - heading and table font sizes
             */
            .manuscript-summary,
            .manuscript-summary * {
              font-size: 12pt !important;
            }

            .manuscript-summary {
              font-family: Arial, sans-serif;
              line-height: 1.25;
              color: #000000 !important;

              -webkit-print-color-adjust: exact !important;
              print-color-adjust: exact !important;
            }

            /*
             * Avoid splitting individual cards, alerts, and table rows
             * across pages where possible.
             */
            .manuscript-summary .card,
            .manuscript-summary .alert,
            .manuscript-summary tr {
              break-inside: avoid;
              page-break-inside: avoid;
            }

            /* Hide interactive elements from the PDF. */
            .manuscript-summary button,
            .manuscript-summary .btn,
            .manuscript-summary .popover,
            .manuscript-summary .tooltip,
            .no-print {
              display: none !important;
            }
          }
          "
        )
      )
    ),
    
    shinyjs::useShinyjs(),
    
    titlePanel(
      "MetaboRx: Quality Assessment and Preprocessing for Metabolomics Data"
    ),
    
    bslib::navset_tab(
      id = "main_steps",
      mod_import_ui("import"),
      mod_correct_ui("correct"),
      mod_visualize_ui("viz"),
      mod_export_ui("export")
    )
  )
}