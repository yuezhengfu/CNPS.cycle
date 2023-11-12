#' ko abundance table related to nitrogen element was extracted
#' @param kegg.entry A profile table containing kegg's entry
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#'
#' N.abundance = Ncyc.abundance(kegg.entry = ko)
#'
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export


Ncyc.abundance <- function(kegg.entry = ko){
    # List of KEGG entry identifiers
    N <- c("K00370","K00371","K00374","K02567","K02568","K00362","K00363","K03385",
           "K15876","K00367","K00372","K00366","K00368","K15864","K04561","K02305",
           "K00376","K02588","K02586","K02591","K10944","K10945","K10946","K10535",
           "K05601","K20932","K20933","K20934","K20935","K01915","K00265","K00266",
           "K00264","K00284","K00260","K15371","K00261","K00262","K01455","K01501",
           "K01725","K00926","K00549","K02575","K15576","K15577","K15578","K15579")

    # Initialize an empty data frame for KEGG entries
    N1 <- c()
    # Populate N1 with data from 'ko' if available
    for (i in 1:length(N)) {
        if (N[i] %in% rownames(ko)) {
            N2 <- ko[N[i],]
        }else{
            N2 <- rep(0,ncol(ko))
        }
        N1 <- rbind(N1,N2)
    }
    rownames(N1) <- N
    N1 <- as.data.frame(t(N1))
    # Calculate N.abundance based on N1
    N.abundance <- data.frame(DNR.NO3.NO2 = (N1$K00370 + N1$K00371 + N1$K00374)/3 +
                                  (N1$K02567 + N1$K02568)/2,
                              DNR.NO2.NH4 = (N1$K00362 + N1$K00363)/2 +
                                  (N1$K03385 + N1$K15876)/2,
                              ANR.NO3.NO2 = N1$K00367 + N1$K00372,
                              ANR.NO2.NH4 = N1$K00366,
                              Denitrification.NO2.NO = N1$K00368 + N1$K15864,
                              Denitrification.NO.N2O = (N1$K04561 + N1$K02305)/2,
                              Denitrification.N2O.N2 = N1$K00376,
                              NF.N2.NH4 = (N1$K02588 + N1$K02586 + N1$K02591)/3,
                              Nitrification.NH4.NH2OH = (N1$K10944 + N1$K10945 + N1$K10946)/3,
                              Nitrification.NH2OH.NO2 = N1$K10535,
                              Nitrification.NH2OH.NH4 = N1$K05601 + N1$K15864,
                              Anammox.NH4.N2H4 = (N1$K20932 + N1$K20933 + N1$K20934)/3,
                              Anammox.N2H4.N2 = N1$K20935,
                              NA.NH4.ON = N1$K01915 + (N1$K00265 + N1$K00266)/2 +
                                  N1$K00264 + N1$K00284,
                              NM.ON.NH4 = N1$K00260 + N1$K15371 + N1$K00261 +
                                  N1$K00262 + N1$K01455 + N1$K01501 + N1$K01725 +
                                  N1$K00926 + N1$K00549,
                              NU.NH4 = N1$K02575,
                              NU.NO2 = (N1$K15576 + N1$K15577 + N1$K15578 + N1$K15579)/4,
                              NU.NO3 = (N1$K15576 + N1$K15577 + N1$K15578 + N1$K15579)/4)
    # Set row names for N.abundance
    rownames(N.abundance) <- rownames(N1)
    # Add a 'Type' column with meaningful names
    N.abundance <- as.data.frame(t(N.abundance))
    N.abundance$Type <- rownames(N.abundance)
    N.abundance <- N.abundance[,c("Type",colnames(N.abundance)[1:(ncol(N.abundance)-1)])]
    N.abundance$Type <- c("Dissimilatory nitrate reduction, nitrate to nitrite",
                          "Dissimilatory nitrate reduction, nitrite to ammonia",
                          "Assimilatory nitrate reduction, nitrate to nitrite",
                          "Assimilatory nitrate reduction, nitrite to ammonia",
                          "Denitrification, nitrite to nitric oxide",
                          "Denitrification, nitric oxide to nitrous oxide",
                          "Denitrification, nitrous oxide to nitrogen",
                          "Nitrogen fixation, nitrogen to ammonia",
                          "Nitrification, ammonia to hydroxylamine",
                          "Nitrification, hydroxylamine to nitrite",
                          "Nitrification, hydroxylamine to ammonia",
                          "Anammox, ammonia to hydrazine",
                          "Anammox, hydrazine to nitrogen",
                          "Nitrogen assimilation, ammonia to organ-N",
                          "Nitrogen mineralization, organ-N to ammonia",
                          "Ammonia uptake",
                          "Nitrite uptake",
                          "Nitrate uptake")

    N.abundance$Type <- factor(N.abundance$Type,levels = rev(N.abundance$Type))
    return(N.abundance)
}
