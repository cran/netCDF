\name{read.netCDF}
\alias{read.netCDF}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Read subsets of a netCDF file}
\description{
Reads a hyperrectangle of values for the specified variables from a file
in netCDF format. 
}
\usage{
read.netCDF(x, name=NULL, id=NULL, start=NULL, count=NULL, byrow=T, attr=T)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{A filename or a \code{netCDF} object as returned by \link{open.netCDF}}
 \item{name}{Variable name(s) to read}
 \item{id}{Variable id(s) to read (specify this or \code{name})}
 \item{start}{upper left corner of hyperrectangle (or list of these if
   reading more than one variable)}
 \item{count}{dimensions of hyperrectangle (or list of these if reading
   more than one variable}}
 \item{byrow}{Transpose the data so that the order of dimensions in R is
 the same as in the file}
 \item{attr}{Return attributes of the variables as well}
}
\details{
  If a filename is given it will be opened before reading and closed
  afterwards. If a \code{netCDF} object is given it will be reopened if
  necessary but will not be closed.
  
  Variables can be described by their name or their variable id. If only
  one variable is being read, or if all variable have the same
  dimensions, \code{start} and \code{count} are vectors whose length is
  the number of dimensions. Otherwise they should be lists of vectors.

  NetCDF files store data with the first dimension changing slowest, the
  opposite to R. The \code{byrow=T} transposes the data so that the
  first dimension changes fastest. This is easier to understand, but
  presumably takes more time and memory.
  
}
\value{
  A matrix or list of matrices containing the data, with attributes if \code{attr=T}
}
\references{Unidata NetCDF website \url{http://www.unidata.ucar.edu/packages/netcdf/}}
\author{Thomas Lumley}
\note{This code uses the deprecated netCDF version 2 compatibility
  library and should be updated.}


\seealso{ \code{\link{open.netCDF}}, \code{\link{close.netCDF}} }

\examples{

}
\keyword{file}%-- one or more ...
