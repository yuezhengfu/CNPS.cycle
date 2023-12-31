% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SC.H2S.Sn.R
\name{SC.H2S.Sn}
\alias{SC.H2S.Sn}
\title{SC.H2S.Sn abundance information was extracted and combined at different species taxonomic levels}
\usage{
SC.H2S.Sn(
  kegg.category = Gene,
  NR.tax = tax,
  profile = abundance,
  group = group
)
}
\arguments{
\item{kegg.category}{kegg grouping file}

\item{NR.tax}{A result of species annotation of gene sets using NR database}

\item{profile}{Results of abundance of gene sets}

\item{group}{Sample grouping file}
}
\description{
SC.H2S.Sn abundance information was extracted and combined at different species taxonomic levels
}
\examples{
data(ko)
data(group)
data(Gene)
data(tax)
data(abundance)

SC.H2S.Sn <- SC.H2S.Sn(kegg.category = Gene,NR.tax = tax,profile = abundance,group = group)

}
\author{
contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
}
