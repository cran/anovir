% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_descriptions.R
\docType{data}
\name{data_blanford}
\alias{data_blanford}
\title{Full data from Blanford et al (2012)}
\format{
An object of class \code{data.frame} with 1320 rows and 9 columns.
}
\source{
Simon Blanford and Matthew Thomas generously provided and allowed the
  release of these data
}
\usage{
data_blanford
}
\description{
Complete data from the publication: Blanford S, Jenkins NE, Read AF, Thomas
MB (2012) Evaluating the lethal and pre-lethal effects of a range of fungi
against adult Anopheles stephensi mosquitoes. Malaria Journal. 11:365
\href{https://doi.org/10.1186/1475-2875-11-365}{doi}
}
\details{
\describe{
  \item{block}{experimental block within experiment (1 - 5)}
  \item{treatment}{experimental treatment}
  \item{replicate cage}{replicate cage within treatment (1 - 4)}
  \item{day}{time post-infection (days)}
  \item{censor}{'1' censored, '0' died }
  \item{d}{an indicator variable; '0' censored, '1' died}
  \item{inf}{'0' uninfected treatement, '1' infected treatment}
  \item{t}{time post-infection (days)}
  \item{fq}{frequency of individuals}
  }
}
\examples{
head(data_blanford)
}
\keyword{datasets}
