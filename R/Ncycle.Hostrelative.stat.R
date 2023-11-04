#' Determine whether a set of genes belongs to a specific biogeochemical cycle
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#' result <- Ncyc.host(Gene)
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com}
#' @export

Ncyc.host <- function(Gene){
  Ncyc <- list(DNR.NO3.NO2 = c("K00370","K00371","K00374","K02567","K02558"),
               DNR.NO2.NH4 = c("K00362","K00363","K03385","K15876"),
               ANR.NO3.NO2 = c("K00367","K00372"),
               ANR.NO2.NH4 = c("K00366"),
               Denitrification.NO2.NO = c("K00368","K15864"),
               Denitrification.NO.N2O  = c("K04561","K02305"),
               Denitrification.N2O.N2 = c("K00376"),
               NF.N2.NH4 = c("K02588","K02586","K02591"),
               Nitrification.NH4.NH2OH = c("K10944","K10945","K10946"),
               Nitrification.NH2OH.NO2  = c("K10535"),
               Nitrification.NH2OH.NH4 = c("K05601","K15864"),
               Anammox.NH4.N2H4 = c("K20932","K20933","K20934"),
               Anammox.N2H4.N2 = c("K20935"),
               NA.NH4.ON = c("K01915","K00265","K00266","K00264","K00284"),
               NM.ON.NH4 = c("K00260","K15371","K00261","K00262","K01455","K01501","K01725","K00926","K00549"),
               NU.NH4 = c("K02575"),
               NU.NO2 = c("K15576","K15577","K15578","K15579"),
               NU.NO3 = c("K15576","K15577","K15578","K15579"))
  temp <- c()
  for (i in 1:18) {
    d <- Ncyc[[i]]
    d <- ifelse(d %in% colnames(Gene),d,NA)
    d <- d[complete.cases(d)]
    a <- ifelse(length(d) == 0,0,1)
    temp <- c(temp,a)
  }
  return(temp)
}
