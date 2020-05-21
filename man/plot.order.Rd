\name{plot.order}
\alias{plot.order}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots ordered datasets
}
\description{
In the process of normalization, statistical depth based or quantile, each element of the sample is ordered in ascending order. This function plots the sample elements in this process. 
}
\usage{
plot.order(ar, x = 1:nrow(ar), tp = "l", xl = range(x), yl = c(min(ar[x, ]), max(ar[x, ])), xlb = "", ylb = "", m = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{ar}{
  Cell files from Affy micrarray experiment: a matrix where the columns are the sample elements, RNA-seq or microarray, and the rows the variables (the genes).
}
  \item{x}{
A vector that is subset of the natural numbers up to the number of variables. It represents the domain, of the X-axis, to use for the zoomed plot. The default is the whole domain.
}
  \item{tp}{
The type of plot that should be drawn. The default is "l" for lines.
}
  \item{xl}{
A vector with the limits of the X-axis. The default is the range of "x".
}
  \item{yl}{
A vector with the limits of the Y-axis. The default is the range of "ar" when restricted to the domain "x".
}
  \item{xlb}{
X-axis label for the plot. The default is empty.
}
  \item{ylb}{
Y-axis label for the plot. The default is empty.
}
  \item{m}{
Main title plot. The default is empty.
}
  \item{\dots}{
Additional arguments to be passed to the plot, such as graphical parameters.
}
}
\details{
This function plots allows to plot a zoom of the domain. The color shows the deepness of the sample elemnt with respect to the sample, using \emph{prof.funct}. The colours go from blue meaning low depth to red meaning high depth, through green and yellow.
}
\value{
}
\references{
Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene expression data. Preprint.
}
\author{
Alicia Nieto-Reyes and Javier Cabrera
}

\seealso{
See Also \code{\link{normalization}} and \code{\link{prof.funct}}
}
\examples{
# Plots the Sialin data set after 6 hours of gestation, in the ascending order used in the process of normalization
plot.order(Sialin6)

# A zoom of previous plot in the interval [20000 , 30000]
plot.order(Sialin6, c(20000:30000))

}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
