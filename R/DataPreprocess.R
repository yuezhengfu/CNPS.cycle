#' Data Preprocessing for kegg.profile.entry, group, kegg.category, NR.taxonomy, profile
#' @param kegg.entry A profile table containing kegg's entry
#' @param group Sample grouping file
#' @param kegg.category kegg grouping file
#' @param NR.tax A result of species annotation of gene sets using NR database
#' @param profile Results of abundance of gene sets
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#' Data.pre(kegg.entry = ko,group = group,kegg.category = Gene,NR.tax = tax,profile = abundance)
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export


Data.pre <- function(kegg.entry = ko,group = group,kegg.category = Gene,NR.tax = tax,profile = abundance) {
  #Remove the last column from 'ko'
  ko <- ko[,-ncol(ko)]
  #Set row names to the first row of 'ko'
  rownames(ko) <- ko[,1]
  # Remove the first column from 'ko'
  ko <- ko[,-1]
  # Convert all columns in 'ko' to numeric
  for (i in 1:ncol(ko)) {
    ko[,i] <- as.numeric(ko[,i])
  }
  # Set column names for 'group'
  colnames(group) <- c("ID","Group")
  # Calculate the number of unique groups and samples
  Group_numb <- length(unique(group[,2]))
  Sample_numb <- length(unique(group[,1]))
  # Remove duplicate rows based on columns 1 and 2 in 'Gene'
  Gene <- Gene[!duplicated(Gene[,1:2]),1:2]
  Gene <- Gene %>%
    group_by(Entry) %>%
    mutate(index = row_number()) %>%
    pivot_wider(names_from = Entry,
                values_from = GeneID) %>%
    select(-index)
  # Select columns 1 and 7 in 'tax' and rename them
  colnames(tax) <- c("V1","V2")
  # Modify 'V2' column values
  tax$V2 <- gsub(".*;k__","",tax$V2)
  tax$V2 <- paste("k__",tax$V2)
  # Set row names for 'abundance'
  rownames(abundance) <- abundance[,1]
  # Remove the first column from 'abundance'
  abundance <- abundance[,-1]
  # Convert all columns in 'abundance' to numeric
  for (i in 1:ncol(abundance)) {
    abundance[,i] <- as.numeric(abundance[,i])
  }
  abundance$V1 <- rownames(abundance)
  abundance <- abundance[,c("V1",colnames(abundance)[1:(ncol(abundance)-1)])]
}

