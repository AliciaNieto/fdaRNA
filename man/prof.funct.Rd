\name{prof.funct}
\alias{prof.funct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Emprirical statistical depth
}
\description{
Computes the sample statistical depth of a given sample using the procedure described in Nieto-Reyes (2011) and Nieto-Reyes and Cabrera (2020).
}
\usage{
prof.funct(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix where the columns are the sample elements and the rows the variables. The matrix can represent functional or multivariate data. It requires at least two variables. When applied to gene expression data, each column is a RNA-seq or microarray and the rows represent the genes. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a list containing the following 
components:
 \item{depth.va}{
 A vector with the depth values of the sample elements.}
 \item{deepest.ele}{
 The identification number (column) of the deepest(s) elements.}
 \item{distant.val}{
 In the process of computing the depth values, these are ordered in pairs. There is a number of pairs equal to the integer part of the sample size divided by two. This vector contains the distance among each two elements of a pair, in decreasing order.}
}
\references{
Nieto-Reyes A. (2011) On the Properties of Functional Depth. In: Ferraty F. (eds) Recent Advances in Functional Data Analysis and Related Topics. Contributions to Statistics. Physica-Verlag HD.

Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene expression data. Preprint.
}
\author{
Alicia Nieto-Reyes and Javier Cabrera
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Applies the "prof.funct" function to the Tissue dataset
p = prof.funct(Tissue)

# Gives the depth values of each of the 41 microarrays in the Tissue dataset
p$depth.val

# Gives the deepest microarray of the Tissue dataset, the one with highest depth value.
p$deepest.ele
}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
