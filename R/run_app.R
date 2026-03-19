#' Launch the ICB Response Explorer Shiny app
#'
#' @param data_dir Optional path to a data directory that overrides the bundled
#'   package data under `inst/app/data`.
#' @param launch.browser Logical; passed to [shiny::runApp()].
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return The return value from [shiny::runApp()].
#' @export
run_app <- function(data_dir = NULL, launch.browser = TRUE, ...) {
  required_pkgs <- c(
    "shiny", "shinyjqui", "tidyverse", "data.table", "DT",
    "scales", "ggpubr", "memoise", "cachem"
  )
  missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "),
      ". Install them before running run_app().",
      call. = FALSE
    )
  }

  app_dir <- system.file("app", package = "ICBResponse")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    dev_app_dir <- normalizePath(file.path(getwd(), "inst", "app"), mustWork = FALSE)
    if (dir.exists(dev_app_dir)) {
      app_dir <- dev_app_dir
    } else {
      stop("Could not locate bundled app directory.", call. = FALSE)
    }
  }

  resolved_data_dir <- data_dir
  if (is.null(resolved_data_dir) || !nzchar(resolved_data_dir)) {
    resolved_data_dir <- getOption("ICBResponse.data_dir", default = "")
  }
  if (!nzchar(resolved_data_dir)) {
    resolved_data_dir <- Sys.getenv("ICBResponse_DATA_DIR", unset = "")
  }
  if (nzchar(resolved_data_dir)) {
    resolved_data_dir <- normalizePath(resolved_data_dir, winslash = "/", mustWork = FALSE)
  }

  if (!nzchar(resolved_data_dir)) {
    app_data_dir <- file.path(app_dir, "data")
    if (dir.exists(app_data_dir)) {
      resolved_data_dir <- normalizePath(app_data_dir, winslash = "/", mustWork = TRUE)
    } else {
      dev_data_dir <- normalizePath(file.path(getwd(), "inst", "app", "data"), winslash = "/", mustWork = FALSE)
      if (dir.exists(dev_data_dir)) {
        resolved_data_dir <- dev_data_dir
      }
    }
  }

  old_data_dir <- getOption("ICBResponse.data_dir")
  on.exit(options(ICBResponse.data_dir = old_data_dir), add = TRUE)
  if (!nzchar(resolved_data_dir)) {
    stop(
      paste(
        "Could not locate the ICBResponse data directory.",
        "Reinstall the package with bundled `inst/app/data`,",
        "or supply `data_dir` explicitly."
      ),
      call. = FALSE
    )
  }

  options(ICBResponse.data_dir = resolved_data_dir)

  shiny::runApp(appDir = app_dir, launch.browser = launch.browser, ...)
}

