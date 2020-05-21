\name{Khan}
\alias{Khan}
\docType{data}
\title{
Four tumors microarray data
}
\description{
The Khan dataset consists of gene expression measurements, obtained using cDNA microarrays, of four types of pediatric small round blue cell tumors: 
the Ewing family of tumors (EW, 23 microarrays), Burkitt lymphomas (BL, 8 microarrays), neuroblastoma (NB, 12 microarrays) and rhabdomyosarcoma (RM, 20 microarrays). There is a total of  63 microarrays and an additional set of 25 microarys for testing purposes. It has 6,567 genes that were filtered to 2,308 genes in Khan et al. (2001). Here it is provided the filtered sample. 
}
\usage{data("Khan")}
\format{
  A data frame with 2308 observations on the following 88 variables.
  \describe{
    \item{\code{TRAIN1.EW}}{a numeric vector}
    \item{\code{TRAIN2.EW}}{a numeric vector}
    \item{\code{TRAIN3.EW}}{a numeric vector}
    \item{\code{TRAIN4.EW}}{a numeric vector}
    \item{\code{TRAIN5.EW}}{a numeric vector}
    \item{\code{TRAIN6.EW}}{a numeric vector}
    \item{\code{TRAIN7.EW}}{a numeric vector}
    \item{\code{TRAIN8.EW}}{a numeric vector}
    \item{\code{TRAIN9.EW}}{a numeric vector}
    \item{\code{TRAIN10.EW}}{a numeric vector}
    \item{\code{TRAIN11.EW}}{a numeric vector}
    \item{\code{TRAIN12.EW}}{a numeric vector}
    \item{\code{TRAIN13.EW}}{a numeric vector}
    \item{\code{TRAIN14.EW}}{a numeric vector}
    \item{\code{TRAIN15.EW}}{a numeric vector}
    \item{\code{TRAIN16.EW}}{a numeric vector}
    \item{\code{TRAIN17.EW}}{a numeric vector}
    \item{\code{TRAIN18.EW}}{a numeric vector}
    \item{\code{TRAIN19.EW}}{a numeric vector}
    \item{\code{TRAIN20.EW}}{a numeric vector}
    \item{\code{TRAIN21.EW}}{a numeric vector}
    \item{\code{TRAIN22.EW}}{a numeric vector}
    \item{\code{TRAIN23.EW}}{a numeric vector}
    \item{\code{TRAIN24.BL}}{a numeric vector}
    \item{\code{TRAIN25.BL}}{a numeric vector}
    \item{\code{TRAIN26.BL}}{a numeric vector}
    \item{\code{TRAIN27.BL}}{a numeric vector}
    \item{\code{TRAIN28.BL}}{a numeric vector}
    \item{\code{TRAIN29.BL}}{a numeric vector}
    \item{\code{TRAIN30.BL}}{a numeric vector}
    \item{\code{TRAIN31.BL}}{a numeric vector}
    \item{\code{TRAIN32.NB}}{a numeric vector}
    \item{\code{TRAIN33.NB}}{a numeric vector}
    \item{\code{TRAIN34.NB}}{a numeric vector}
    \item{\code{TRAIN35.NB}}{a numeric vector}
    \item{\code{TRAIN36.NB}}{a numeric vector}
    \item{\code{TRAIN37.NB}}{a numeric vector}
    \item{\code{TRAIN38.NB}}{a numeric vector}
    \item{\code{TRAIN39.NB}}{a numeric vector}
    \item{\code{TRAIN40.NB}}{a numeric vector}
    \item{\code{TRAIN41.NB}}{a numeric vector}
    \item{\code{TRAIN42.NB}}{a numeric vector}
    \item{\code{TRAIN43.NB}}{a numeric vector}
    \item{\code{TRAIN44.RM}}{a numeric vector}
    \item{\code{TRAIN45.RM}}{a numeric vector}
    \item{\code{TRAIN46.RM}}{a numeric vector}
    \item{\code{TRAIN47.RM}}{a numeric vector}
    \item{\code{TRAIN48.RM}}{a numeric vector}
    \item{\code{TRAIN49.RM}}{a numeric vector}
    \item{\code{TRAIN50.RM}}{a numeric vector}
    \item{\code{TRAIN51.RM}}{a numeric vector}
    \item{\code{TRAIN52.RM}}{a numeric vector}
    \item{\code{TRAIN53.RM}}{a numeric vector}
    \item{\code{TRAIN54.RM}}{a numeric vector}
    \item{\code{TRAIN55.RM}}{a numeric vector}
    \item{\code{TRAIN56.RM}}{a numeric vector}
    \item{\code{TRAIN57.RM}}{a numeric vector}
    \item{\code{TRAIN58.RM}}{a numeric vector}
    \item{\code{TRAIN59.RM}}{a numeric vector}
    \item{\code{TRAIN60.RM}}{a numeric vector}
    \item{\code{TRAIN61.RM}}{a numeric vector}
    \item{\code{TRAIN62.RM}}{a numeric vector}
    \item{\code{TRAIN63.RM}}{a numeric vector}
    \item{\code{TEST64.NB}}{a numeric vector}
    \item{\code{TEST65.EW}}{a numeric vector}
    \item{\code{TEST66.NA}}{a numeric vector}
    \item{\code{TEST67.RM}}{a numeric vector}
    \item{\code{TEST68.NA}}{a numeric vector}
    \item{\code{TEST69.EW}}{a numeric vector}
    \item{\code{TEST70.BL}}{a numeric vector}
    \item{\code{TEST71.NB}}{a numeric vector}
    \item{\code{TEST72.NA}}{a numeric vector}
    \item{\code{TEST73.RM}}{a numeric vector}
    \item{\code{TEST74.NA}}{a numeric vector}
    \item{\code{TEST75.EW}}{a numeric vector}
    \item{\code{TEST76.NA}}{a numeric vector}
    \item{\code{TEST77.NB}}{a numeric vector}
    \item{\code{TEST78.BL}}{a numeric vector}
    \item{\code{TEST79.NB}}{a numeric vector}
    \item{\code{TEST80.RM}}{a numeric vector}
    \item{\code{TEST81.BL}}{a numeric vector}
    \item{\code{TEST82.EW}}{a numeric vector}
    \item{\code{TEST83.EW}}{a numeric vector}
    \item{\code{TEST84.EW}}{a numeric vector}
    \item{\code{TEST85.RM}}{a numeric vector}
    \item{\code{TEST86.NB}}{a numeric vector}
    \item{\code{TEST87.RM}}{a numeric vector}
    \item{\code{TEST88.NB}}{a numeric vector}
  }
}
\details{
Variables:

Col. 1-23:  	TRAIN1.EW through TRAIN23.EW. Training set with tumor type EW.

Col. 24-31: 	TRAIN24.BL through TRAIN31.BL. Training set with tumor type EW.

Col. 32-43: 	TRAIN32.NB through TRAIN43.NB. . Training set with tumor type EW.

Col. 34-63:	  TRAIN44.RM through TRAIN63.RM. Training set with tumor type EW.

Col. 64-88:   EST64.NB through TEST88.NB. Testing set.
}
\source{
Own
}
\references{
Amaratunga D, Cabrera J, Shkedy Z. Exploration and analysis of DNA microarray and other high dimensional data. J. Wiley & Sons, 2014. 

Khan J, Wei J, Ringner M, Saal L, Ladanyi M, Westermann F, Berthold F, Schwab M, Antonescu C, Peterson C, and Meltzer P. Classification and diagnostic prediction of cancers using gene expression profiling and artificial neural networks. Nature Medicine, v.7, pp.673-679, 2001.

Nieto-Reyes A, Cabrera J. Statistical depth based normalization and outlier detection of gene expression data. Preprint.
}
\examples{
K=Khan[,1:63]
names(K)
outlier_detection(K,MS=T)
}
\keyword{datasets}
