#' Determine whether a set of genes belongs to a specific biogeochemical cycle
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#'
#' result <- Ccyc.host(Gene)
#'
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export

Ccyc.host <- function(Gene){
  Ccyc <- list(ACF = c("K00855","K01061","K01062"),
               ACH4O = c("K08684"),
               AR = c("K02256","K02262","K02274","K02276"),
               AnCF = c("K00174","K00175","K00244","K01648","K00194","K00197"),
               COo = c("K03518","K03519","K03520"),
               Fer = c("K00016"),
               Meth = c("K00400","K00401"))
  temp <- c()
  for (i in 1:7) {
    d <- Ccyc[[i]]
    d <- ifelse(d %in% colnames(Gene),d,NA)
    d <- d[complete.cases(d)]
    a <- ifelse(length(d) == 0,0,1)
    temp <- c(temp,a)
  }
  return(temp)
}
