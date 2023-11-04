#' ko abundance table related to sulfur element was extracted
#' @param kegg.entry A profile table containing kegg's entry
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#' S.abundance = Scyc.abundance(kegg.entry = ko)
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com}
#' @export



Scyc.abundance <- function(kegg.entry = ko){
  # List of KEGG entry identifiers
    S <- c("K00956","K00957","K00955","K00860","K00390","K00380","K00381","K00392",
           "K00958","K00988","K00394","K00395","K11180","K11181","K17222","K17223",
           "K17224","K17225","K17226","K17227","K17229","K17230","K17218","K08352",
           "K08354","K04091","K00299","K16968","K16969","K15762","K15765","K03119",
           "K00456","K17217","K02045","K02046","K02047","K02048","K15551","K15552",
           "K10831","K15553","K15554","K15555","K01739","K10764","K01738","K17228")
  # Initialize an empty data frame for KEGG entries
    S1 <- c()
  # Populate S1 with data from 'ko' if available
    for (i in 1:length(S)) {
        if (S[i] %in% rownames(ko)) {
            S2 <- ko[S[i],]
        }else{
            S2 <- rep(0,ncol(ko))
        }
        S1 <- rbind(S1,S2)
    }
    rownames(S1) <- S
    S1 <- as.data.frame(t(S1))

  # Calculate S.abundance based on S1
    S.abundance <- data.frame(ASR.SO4.APS = (S1$K00956 + S1$K00957)/2 + S1$K00955,
                              ASR.APS.PAPS = S1$K00955 + S1$K00860,
                              ASR.PAPS.SO3 = S1$K00390,
                              ASR.SO3.H2S = (S1$K00380 + S1$K00381)/2 + S1$K00392,
                              DSR.SO4.APS = S1$K00958 + S1$K00988,
                              DSR.APS.SO3 = (S1$K00394 + S1$K00395)/2,
                              DSR.SO3.H2S = (S1$K11180 + S1$K11181)/2,
                              SOX = (S1$K17222 + S1$K17223)/2 + S1$K17224 +
                                  S1$K17225 + (S1$K17226 + S1$K17227)/2,
                              SC.H2S.S = (S1$K17229 + S1$K17230)/2,
                              SC.H2S.Sn = S1$K17218,
                              SC.T.H2S = (S1$K08352 + S1$K08354)/2,
                              SM = (S1$K11180 + S1$K11181)/2 + S1$K04091 + S1$K00299 +
                                  S1$K17228 + (S1$K16968 + S1$K16969)/2 +
                                  (S1$K15762 + S1$K15765)/2 + S1$K03119 +
                                  S1$K00456 + S1$K17217,
                              SU.SO4 = (S1$K02045 + S1$K02046 + S1$K02047 +
                                            S1$K02048)/4,
                              SU.SO3 = (S1$K15551 + S1$K15552 + S1$K10831)/3 +
                                  (S1$K15554 +S1$K15555 + S1$K15553)/3,
                              SA = S1$K01739 + S1$K10764 + S1$K01738)

    # Set row names for S.abundance
    rownames(S.abundance) <- rownames(S1)
    # Add a 'Type' column with meaningful names
    S.abundance <- as.data.frame(t(S.abundance))
    S.abundance$Type <- rownames(S.abundance)
    S.abundance <- S.abundance[,c("Type",colnames(S.abundance)[1:(ncol(S.abundance)-1)])]
    S.abundance$Type <- c("Assimilatory sulfate reduction, sulfate to APS",
                          "Assimilatory sulfate reduction, APS to PAPS",
                          "Assimilatory sulfate reduction, PAPS to sulfite",
                          "Assimilatory sulfate reduction, sulfite to sulfide",
                          "Dissimilatory sulfate reduction and oxidation, sulfate to APS",
                          "Dissimilatory sulfate reduction and oxidation, APS to sulfite",
                          "Dissimilatory sulfate reduction and oxidation, sulfite to sulfide",
                          "SOX system",
                          "Sulfide cycling, sulfide to sulfur",
                          "Sulfide cycling, sulfide to (sulfide)n",
                          "Sulfide cycling, thisulfate to sulfide",
                          "Sulfur mineralization",
                          "Sulfate uptake",
                          "Sulfite uptake",
                          "Sulfur assimilation")
    # Reorder the 'Type' column
    S.abundance$Type <- factor(S.abundance$Type,levels = rev(S.abundance$Type))
    return(S.abundance)
}
