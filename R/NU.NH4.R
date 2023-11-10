#' NU.NH4 abundance information was extracted and combined at different species taxonomic levels
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
#' NU.NH4 <- NU.NH4(kegg.category = Gene,NR.tax = tax,profile = abundance,group = group)
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export

NU.NH4 <- function(kegg.category = Gene,NR.tax = tax,profile = abundance,group = group){
    A <- as.data.frame(Gene$K02575)
    A[A == ""] <- NA
    A <- as.data.frame(na.omit(A))
    colnames(A) <- c("V1")
    A <- merge(A,tax,by = "V1")
    C <- merge(A,abundance,by = "V1")

    B <- as.data.frame(C$V2)
    C1 <- as.data.frame(t(C[,3:ncol(C)]))
    C1$ID <- rownames(C1)
    C1 <- merge(C1,group)
    C2 <- aggregate(C1[,2:(ncol(C1)-1)],list(C1$Group),mean)
    rownames(C2) <- C2$Group.1
    C2 <- as.data.frame(C2[,-1])
    if (ncol(C2) == 1) {
        rownames(C2) <- unique(group$Group)
    }
    if (ncol(C2) > 1) {
        C2 <- C2/rowSums(C2)
    }
    if (ncol(C2) == 1) {
        C2 <- C2/C2
    }
    B <- cbind(B,t(C2))
    B <- separate(B,col = `C$V2`,sep = ";",into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

    B$Phylum <- factor(B$Phylum)
    group_id <- colnames(B)[8:ncol(B)]
    B1 <- c()
    for (i in group_id) {
        aa <- tapply(B[,i], B$Phylum, sum)
        B1 <- cbind(B1,aa)
    }
    colnames(B1) <- group_id
    B1 <- as.data.frame(B1)
    B1$Sum <- apply(B1,1,sum)
    B1 <- B1[order(B1[,ncol(B1)],decreasing = TRUE),]
    B1 <- subset(B1,select = -Sum)
    if("p__Unclassified" %in% rownames(B1)){
        B1 <- t(B1)
        B1 <- subset(B1,select = -p__Unclassified)
        B1 <- t(B1)
    }
    rownames(B1) <- sub("^...","",rownames(B1))
    Phylum <- B1

    B$Class <- factor(B$Class)
    group_id <- colnames(B)[8:ncol(B)]
    B1 <- c()
    for (i in group_id) {
        aa <- tapply(B[,i], B$Class, sum)
        B1 <- cbind(B1,aa)
    }
    colnames(B1) <- group_id
    B1 <- as.data.frame(B1)
    B1$Sum <- apply(B1,1,sum)
    B1 <- B1[order(B1[,ncol(B1)],decreasing = TRUE),]
    B1 <- subset(B1,select = -Sum)
    if("c__Unclassified" %in% rownames(B1)){
        B1 <- t(B1)
        B1 <- subset(B1,select = -c__Unclassified)
        B1 <- t(B1)
    }
    rownames(B1) <- sub("^...","",rownames(B1))
    Class <- B1

    B$Order <- factor(B$Order)
    group_id <- colnames(B)[8:ncol(B)]
    B1 <- c()
    for (i in group_id) {
        aa <- tapply(B[,i], B$Order, sum)
        B1 <- cbind(B1,aa)
    }
    colnames(B1) <- group_id
    B1 <- as.data.frame(B1)
    B1$Sum <- apply(B1,1,sum)
    B1 <- B1[order(B1[,ncol(B1)],decreasing = TRUE),]
    B1 <- subset(B1,select = -Sum)
    if("o__Unclassified" %in% rownames(B1)){
        B1 <- t(B1)
        B1 <- subset(B1,select = -o__Unclassified)
        B1 <- t(B1)
    }
    rownames(B1) <- sub("^...","",rownames(B1))
    Order <- B1

    B$Family <- factor(B$Family)
    group_id <- colnames(B)[8:ncol(B)]
    B1 <- c()
    for (i in group_id) {
        aa <- tapply(B[,i], B$Family, sum)
        B1 <- cbind(B1,aa)
    }
    colnames(B1) <- group_id
    B1 <- as.data.frame(B1)
    B1$Sum <- apply(B1,1,sum)
    B1 <- B1[order(B1[,ncol(B1)],decreasing = TRUE),]
    B1 <- subset(B1,select = -Sum)
    if("f__Unclassified" %in% rownames(B1)){
        B1 <- t(B1)
        B1 <- subset(B1,select = -f__Unclassified)
        B1 <- t(B1)
    }
    rownames(B1) <- sub("^...","",rownames(B1))
    Family <- B1

    B$Genus <- factor(B$Genus)
    group_id <- colnames(B)[8:ncol(B)]
    B1 <- c()
    for (i in group_id) {
        aa <- tapply(B[,i], B$Genus, sum)
        B1 <- cbind(B1,aa)
    }
    colnames(B1) <- group_id
    B1 <- as.data.frame(B1)
    B1$Sum <- apply(B1,1,sum)
    B1 <- B1[order(B1[,ncol(B1)],decreasing = TRUE),]
    B1 <- subset(B1,select = -Sum)
    if("g__Unclassified" %in% rownames(B1)){
        B1 <- t(B1)
        B1 <- subset(B1,select = -g__Unclassified)
        B1 <- t(B1)
    }
    rownames(B1) <- sub("^...","",rownames(B1))
    Genus <- B1

    B$Species <- factor(B$Species)
    group_id <- colnames(B)[8:ncol(B)]
    B1 <- c()
    for (i in group_id) {
        aa <- tapply(B[,i], B$Species, sum)
        B1 <- cbind(B1,aa)
    }
    colnames(B1) <- group_id
    B1 <- as.data.frame(B1)
    B1$Sum <- apply(B1,1,sum)
    B1 <- B1[order(B1[,ncol(B1)],decreasing = TRUE),]
    B1 <- subset(B1,select = -Sum)
    if("s__Unclassified" %in% rownames(B1)){
        B1 <- t(B1)
        B1 <- subset(B1,select = -s__Unclassified)
        B1 <- t(B1)
    }
    rownames(B1) <- sub("^...","",rownames(B1))
    Species <- B1
    result <- list(Phylum,Class,Order,Family,Genus,Species)
    return(result)
}
