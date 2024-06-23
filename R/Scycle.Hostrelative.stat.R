#' Determine whether a set of genes belongs to a specific biogeochemical cycle
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#'
#' result <- Scyc.host(Gene)
#'
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export

Scyc.host <- function(Gene){
  Scyc <- list(ASR.SO4.APS = c("K00955","K00956","K00957"),
               ASR.APS.PAPS = c("K00955","K00860"),
               ASR.PAPS.SO3 = c("K00390"),
               ASR.SO3.H2S = c("K00380","K00381","K00392"),
               DSR.SO4.APS = c("K00958","K00988"),
               DSR.APS.SO3 = c("K00394","K00395"),
               DSR.SO3.H2S = c("K11180","K11181"),
               SOX = c("K17722","K17723","K17724","K17725","K17726","K17727"),
               SC.H2S.S = c("K17729","K17730"),
               SC.H2S.Sn = c("K17218"),
               SC.T.H2S = c("K08352","K08354"),
               SM = c("K11180","K11181","K04091","K00299","K17728","K16968","K16969","K15762",
                      "K15765","K03119","K00456","K17217"),
               SU.SO4 = c("K02045","K02046","K02047","K02048"),
               SU.SO3 = c("K10831","K15551","K15552","K15553","K15554","K15555"),
               SA = c("K01739","K10764","K01738"))
  temp <- c()
  for (i in 1:15) {
    d <- Scyc[[i]]
    d <- ifelse(d %in% colnames(Gene),d,NA)
    d <- d[complete.cases(d)]
    a <- ifelse(length(d) == 0,0,1)
    temp <- c(temp,a)
  }
  return(temp)
}
