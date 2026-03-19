% Generated manually for CRAN checks
\name{run_app}
\alias{run_app}
\title{Launch the ICB Response Explorer Shiny app}
\usage{
run_app(data_dir = NULL, launch.browser = TRUE, ...)
}
\arguments{
\item{data_dir}{Optional path to a data directory that overrides the bundled
package data under \code{inst/app/data}.}
\item{launch.browser}{Logical; passed to \code{shiny::runApp()}.}
\item{...}{Additional arguments passed to \code{shiny::runApp()}.}
}
\value{The return value from \code{shiny::runApp()}.}
\description{Launches the bundled Shiny application. By default it uses the
bundled package data under \code{inst/app/data}. Use \code{data_dir}, the
\code{ICBResponse.data_dir} option, or \code{ICBResponse_DATA_DIR} environment
variable only to override that default.}
\examples{
\dontrun{
ICBResponse::run_app()
}
}

