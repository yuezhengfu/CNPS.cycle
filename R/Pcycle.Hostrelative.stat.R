#' Determine whether a set of genes belongs to a specific biogeochemical cycle
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#' result <- Pcyc.host(Gene)
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export

Pcyc.host <- function(Gene){
  Pcyc <- list(RB = c("K07636","K07657","K07658"),
               X3 = c("K07768","K07776"))
  temp <- c()
  for (i in 1:2) {
    d <- Pcyc[[i]]
    d <- ifelse(d %in% colnames(Gene),d,NA)
    d <- d[complete.cases(d)]
    a <- ifelse(length(d) == 0,0,1)
    temp <- c(temp,a)
  }
  return(temp)
}
