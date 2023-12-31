% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fold change.R
\name{fold.change}
\alias{fold.change}
\title{Find the fold change value and visualize it using a heat map}
\usage{
fold.change(ARG_type, group)
}
\arguments{
\item{ARG_type}{C.abundance, N.abundance, P.abundance, S.abundance}

\item{group}{Sample grouping file}
}
\description{
Find the fold change value and visualize it using a heat map
}
\examples{
data(ko)
data(group)
data(Gene)
data(tax)
data(abundance)

result <- fold.change(C.abundance,group)

}
\author{
contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
}
