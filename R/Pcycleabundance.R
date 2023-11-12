#' ko abundance table related to phosphorus element was extracted
#' @param kegg.entry A profile table containing kegg's entry
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#'
#' P.abundance = Pcyc.abundance(kegg.entry = ko)
#'
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export


Pcyc.abundance <- function(kegg.entry = ko){
  # List of KEGG entry identifiers
    P <- c("K07636","K07657","K07658","K07768","K07776")
  # Initialize an empty data frame for KEGG entries
    P1 <- c()
    for (i in 1:length(P)) {
        if (P[i] %in% rownames(ko)) {
            P2 <- ko[P[i],]
        }else{
            P2 <- rep(0,ncol(ko))
        }
        P1 <- rbind(P1,P2)
    }
    rownames(P1) <- P
    P1 <- as.data.frame(t(P1))
  # Calculate P.abundance based on P1
    P.abundance <- data.frame(RB = P1$K07636 + (P1$K07657 + P1$K07658)/2,
                              X3 = P1$K07768 + P1$K07776)
  # Set row names for P.abundance
    rownames(P.abundance) <- rownames(P1)
  # Add a 'Type' column with meaningful names
    P.abundance <- as.data.frame(t(P.abundance))
    P.abundance$Type <- rownames(P.abundance)
    P.abundance <- P.abundance[,c("Type",colnames(P.abundance)[1:(ncol(P.abundance)-1)])]
    P.abundance$Type <- c("PhoR-PhoB system","SenX3-RegX3 system")
  # Reorder the 'Type' column
    P.abundance$Type <- factor(P.abundance$Type,levels = rev(P.abundance$Type))
    return(P.abundance)
}
