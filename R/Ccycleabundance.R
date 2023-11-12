#' ko abundance table related to carbon element was extracted
#' @param kegg.entry A profile table containing kegg's entry
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#'
#' C.abundance = Ccyc.abundance(kegg.entry = ko)
#'
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export


Ccyc.abundance <- function(kegg.entry = ko) {
  # List of KEGG entry identifiers
  C <- c("K00855","K01061","K01602","K08684","K02256","K02262","K02274","K02276",
         "K00174","K00175","K00244","K01648","K00194","K00197","K03518",
         "K03519","K03520","K00016","K00400","K00401")
  # Initialize an empty data frame for KEGG entries
  C1 <- c()
  # Populate C1 with data from 'ko' if available
  for (i in 1:length(C)) {
    if (C[i] %in% rownames(ko)) {
      C2 <- ko[C[i],]
    }else{
      C2 <- rep(0,ncol(ko))
    }
    C1 <- rbind(C1,C2)
  }
  rownames(C1) <- C
  C1 <- as.data.frame(t(C1))
  # Calculate C.abundance based on C1
  C.abundance <- data.frame(ACF = C1$K00855 + (C1$K01061 + C1$K01602)/2,
                            ACH4O = C1$K08684,
                            AR = (C1$K02256 + C1$K02262)/2 + (C1$K02274 + C1$K02276)/2,
                            AnCF = (C1$K00174 + C1$K00175)/2 + C1$K00244 +
                              C1$K01648 + (C1$K00194 + C1$K00197)/2,
                            COo = C1$K03518 + (C1$K03519 + C1$K03520)/2,
                            Fer = C1$K00016,
                            Meth = (C1$K00400 + C1$K00401)/2)
  # Set row names for C.abundance
  rownames(C.abundance) <- rownames(C1)
  C.abundance <- as.data.frame(t(C.abundance))
  # Add a 'Type' column with meaningful names
  C.abundance$Type <- rownames(C.abundance)
  C.abundance <- C.abundance[,c("Type",colnames(C.abundance)[1:(ncol(C.abundance)-1)])]

  C.abundance$Type <- c("Aerobic C fixation","Aerobic CH4 oxidation",
                        "Aerobic respiration","Anaerobic C fixation",
                        "CO oxidation","Fermentation","Methanogenesis")
  # Convert 'Type' to a factor with custom order
  C.abundance$Type <- factor(C.abundance$Type,levels = rev(C.abundance$Type))
  return(C.abundance)
}

