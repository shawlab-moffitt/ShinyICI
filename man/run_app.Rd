% Generated manually for CRAN checks
\name{run_app}
\alias{run_app}
\title{Launch the ICB Response Explorer Shiny app}
\usage{
run_app(data_dir = NULL, launch.browser = TRUE, ...)
}
\arguments{
\item{data_dir}{Optional path to the external data directory containing the
files required by the Shiny app.}
\item{launch.browser}{Logical; passed to \code{shiny::runApp()}.}
\item{...}{Additional arguments passed to \code{shiny::runApp()}.}
}
\value{The return value from \code{shiny::runApp()}.}
\description{Launches the bundled Shiny application. For CRAN-sized builds, the
large app data should be stored outside the package and supplied with
\code{data_dir} or the \code{ICBResponse.data_dir} option / \code{ICBResponse_DATA_DIR}
environment variable.}
\examples{
\dontrun{
ICBResponse::run_app(data_dir = "path/to/ICBResponse-data")
}
}
