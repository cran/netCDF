\name{open.netCDF}
\alias{open.netCDF}
\alias{is.open.netCDF}
\alias{close.netCDF}
\alias{close}
\alias{dim.netCDF}
\alias{names.netCDF}
\alias{summary.netCDF}
\alias{print.summary.netCDF}
\alias{print.netCDF}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Open a NetCDF data file}
\description{
These functions open, close, and read metadata from a netCDF data
file. For \code{open.netCDF} the whole function name must be used; the
other functions are methods for existing generic functions.
}
\usage{
open.netCDF(filename, verbose=F)
is.open.netCDF(object)
close.netCDF(object)
dim.netCDF(object)
names.netCDF(object)
summary.netCDF(object)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{filename}{name of netCDF file}
 \item{verbose}{Should netCDF error messages be printed?}
 \item{x}{see below}
}
\details{
NetCDF is a portable data format used by the National Center for
Atmosphereric Research and many others.  Its advantages are that it
allows metadata (units, dimensions,etc) to be stored and that arbitrary
hyperrectangle subsets of the data can be read efficiently.

The \code{open.netCDF} function opens a netCDF file, which remains open
until it is explicitly closed.  Data can be read with \code{\link{read.netCDF}}.

}
\value{
  \code{open.netCDF} returns an object of class \code{netCDF},
  \code{dim()} returns the dimensions, \code{names} returns the variable
  names, \code{summary} returns all the attributes. 
}
\references{Unidata netCDF website: \url{http://www.unidata.ucar.edu/packages/netcdf/}
\author{Thomas Lumley}
\note{This code uses the deprecated netCDF version 2 compatibility
  library and should be updated.}


\seealso{\code{\link{read.netCDF}}}

\examples{
}
\keyword{file}

