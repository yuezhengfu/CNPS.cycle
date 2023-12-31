% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCoA.R
\name{pcoa.arg}
\alias{pcoa.arg}
\title{Principal coordinate analysis (PCoA) and multivariate data analysis}
\usage{
pcoa.arg(ARG_sub, group)
}
\arguments{
\item{ARG_sub}{C.abundance[,2:ncol(C.abundance)], N.abundance[,2:ncol(N.abundance)], P.abundance[,2:ncol(P.abundance)], S.abundance[,2:ncol(S.abundance)]}

\item{group}{Sample grouping file}
}
\description{
Principal coordinate analysis (PCoA) and multivariate data analysis
}
\examples{
data(ko)
data(group)
data(Gene)
data(tax)
data(abundance)

result <- pcoa.arg(C.abundance[,2:ncol(C.abundance)],group)

}
\author{
contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
}
