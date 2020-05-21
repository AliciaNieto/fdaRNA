\name{outlier.detection}
\alias{outlier.detection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Outliers Detection Procedure for Microarrays
}
\description{
This function performs an anlytical outlier detection procedure, for a sample of microarray cell files from a microarray experiment with Affymetric chips, to find microarrays that are outliers as a whole. However, the function can be generally used for a sample of multivariate data, dimension larger than one.
}
\usage{
outlier.detection(ar, nor = T, c.T = F, alpha = 0.01, ap = 1, mn = 1:ncol(ar), t.m = "random", MS = F, SM = F, clases = 1, Cla = list(1:ncol(ar)), ce = 1, ms = "Multidimensional Scaling", sm = "Spectral Map", xlb = " ", ylb = " ", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ar}{
Cell files from Affy micrarray experiment: a matrix where the columns are the sample elements and the rows the variables. When applied to gene expression data, each column is a RNA-seq or microarray and the rows represent the genes. It requires at least two variables.
}
  \item{nor}{
If TRUE, the statistical functional depth based normalization procedure \emph{normalization} is performed to the data prior to the outlier detection technique. Set to FALSE when applied to non gene expression data. The default is TRUE.
}
  \item{c.T}{
This is the cut-off value that determines whether an element of the sample is considered or not as potential outlier. It is not required to provide a value for "c.T", the defalut. See the Details Section. 

If c.T is a real number, the procedure checked whether it is smaller than the ratio of the distance between the elements of a pair and the median distance of the pairs. In that case, it considers an element of the pair as a potential outlier.
  
If no value is provided, "c.T" is obtained by drawing a multivariate normal distribution with the dimension, sample size and covariance structure of "ar" and estimating the quantile such that the "alpha" percent of the elements of the drawn sample are considered as outliers. "c.T" is the median of the quantiles after repating this procedure "ap" times.
}
  \item{alpha}{
If no value of "c.T" is provided, this is the percentage used to estimate it. Otherwise, it is not used. The default is .01.
}
  \item{ap}{
If no value of "c.T" is provided, this is the number of times the procedure is run to estimate it. Otherwis, it is not used. The default is 1.
}
  \item{mn}{
If a subsample of the dataset is aimed to study (while using the whole dataset to do the normalization and/or estimating "c.T""), this is the vector of the columns we aim to study. Thus, if "nor" is TRUE, the default, the normalization procedure is applied to the whole dataset. Anagously, if "c.T" is not provided, it is estimated using the whole dataset.
The default is to study the whole sample.  
}
  \item{t.m}{
A character string specifying how the function "rank{base}" treats ties. The "random" method puts these in random order whereas the default, "average", replaces them by their mean, and "max" and "min" replaces them by their maximum and minimum respectively. The default is "random".
}
  \item{MS}{
If TRUE, it plots the result of performing the multidimensional scaling graphical procedure. The default is FALSE. 
}
  \item{SM}{
If TRUE, it plots the result of performing the spectral map graphical procedure. The default is FALSE. 
}
  \item{clases}{
Number of colors to  be used in the  multidimensional scaling and/or spectral map plots. The default is 1. It is only used if "MS" or/and "SM" is/are set to TRUE.
}
  \item{Cla}{
A list with a number "clases" of elements. Each element is a vector with a subset of columns of "mn". Each element of "mn" has to be in one and only one of the elements of the list. The samples corresponding to the first element of the list are plotted in black, the ones corresponding to the seconde element in red, to the third in green, to the fourth in blue, ... It is only used if "MS" or/and "SM" is/are set to TRUE and "clases" is larger than 1.
}
  \item{ce}{
Number indicating the amount by which the plotted numbers should be scaled relative to the default. The default is 1. It is only used if "MS" or/and "SM" is/are set to TRUE.
}
  \item{ms}{
Main title plot in the multidimensional scaling plot. The default is 'Multidimensional Scaling'. It is only used if "MS"  is set to TRUE.
}
  \item{sm}{
Main title plot in the spectral map plot. The default is 'Spectral Map'. It is only used if "SM" is set to TRUE.
}
  \item{xlb}{
X-axis label for the  multidimensional scaling and/or spectral map plots. The default is empty. It is only used if "MS" or/and "SM" is/are set to TRUE.
}
  \item{ylb}{
Y-axis label for the  multidimensional scaling and/or spectral map plots. The default is empty. It is only used if "MS" or/and "SM" is/are set to TRUE.
}
  \item{\dots}{
Additional arguments to be passed to the multidimensional scaling and/or spectral map plots, such as graphical parameters. It is only used if "MS" or/and "SM" is/are set to TRUE.
}
}
\details{
This procedure does not intend to detect outlying parts in the microarray but whether the microarray is an outlier element with respect to the sample. The anlytical outlier detection procedure that performs this function is based on statistical functional depth.
The procedure first applies the statistical depth to rank the sample elements, using \emph{prof.funct}, and it obtains pairs of sample elements with the same rank. A pair is considered to contain a potential outlier if the distance between its elements is larger than a value "c.T" times the median distance of the pairs. When a pair is a potential outlier, teh function flags as outlier the element fardest away from the statistical deepest microarray. 

Toguether with the analytical procedure, this function allows to perform a Multidimensional Scaling and Spectral Maps plots to compare the analytical results with these graphical methodologies.
}
\value{
The function returns a list containing the following: 
components:
 \item{Distance.values}{A matrix with three rows. The number of columns is the integer part of the sample size divided by two. Each column corresponds to a pair of sample elements. The first two rows are the identification of the element in the sample and the third row is the $L_2$ (euclidean) distance between the elements in the pair. The order of the columns follows the statistical depth as in prof.funct. The first column is formed by the two least deep elements of the sample, the second column by the third and fourth least deep elements of the sample and so on. If the sample size is even, the last column represents the two deepest elements. It it is odd, it represents the elements with second and third highest depth.}
 \item{benchmark.outliers}{The value that, if exceed by a value on the third row of the matrix "Distance.values", an element of the corresponding sample pair is considered as an outlier.}
 \item{$id.out}{If integer(0), no outlier is detected. Otherwise, it contains the identification number of the sample size elements that are considered as outliers.}
 \item{constant.Tukey}{It is the value c.T., it is given. Otherwise it is the corresponding value estimated by means of drawing multivariate normal distribution with the dimension, sample size and covariance structure of "ar".}
}
\references{
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
See Also \code{\link{normalization}} and \code{\link{prof.funct}}
}
\examples{
# Applies the analytical outlier detection procedure on the treated elements of the Airway data, after applying the statistical depth based normalization to the whole dataset. "c.T" is estimated using the whole Airway dataset.
outlier.detection(Airway,mn=c(2, 4, 6 ,8))

# Applies the analytical outlier detection procedure on the treated elements of the Airway data, after applying the statistical depth based normalization to those elements only. "c.T" is estimated using only the treated elements of the Airway dataset.
outlier.detection(Airway[,c(2, 4, 6 ,8)])

# Applies the analytical outlier detection procedure to the Airway data. "c.T" is estimated as the median of the result of running the estimating procedure 100 times.
outlier.detection(Airway, ap=100)

# Applies the analytical outlier detection procedure to the Airway data. "c.T" is estimated as the median of the result of running the estimating procedure 10 times. The Multidimensional Scaling and Spectral Map plots are displayed. The untreated data, control group, is colored in black and the treated data in red.
outlier.detection(Airway, ap=10, clases=2, Cla=list(c(1,3,5,7),c(2, 4, 6 ,8)),MS=T,SM=T)

}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
