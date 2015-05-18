\name{Pheno}
\alias{Pheno}
\docType{data}
\title{
An example of Pheno objects
}
\description{
An object consisting of 100 samples and three traits, "Height"", "Weight"", and "Length". 
}
\details{
Five QTLs are dispersed across the genomes. The QTLs are not shared among traits. The QTL effects were drawn from a normal distribution with mean = 0 and SD = 4. The heritabilities of the traits are 0.8 (Height), 0.5 (Weight), and 0.7 (Length). The "Weight"" record of the fifth sample is missing (NA).
}
\examples{
data(sampledata)
dim(Pheno)#100 samples and 3 traits
any(is.na(Pheno))#TRUE. Pheno includes missing records
which(is.na(Pheno$Weight))#5. The fifth sample lacks the Weight record.
}
\keyword{datasets}
