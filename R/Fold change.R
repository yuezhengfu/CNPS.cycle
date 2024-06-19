#' Find the fold change value and visualize it using a heat map
#' @param ARG_type C.abundance, N.abundance, P.abundance, S.abundance
#' @param group Sample grouping file
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#' result <- fold.change(C.abundance,group)
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export

fold.change <- function(ARG_type,group){
    ARG_type2 <- ARG_type[which(rowSums(ARG_type[,2:ncol(ARG_type)])>0),]
    ARG_type2 <- melt(ARG_type2)
    colnames(ARG_type2) <- c("Type","variable","value")
    colnames(group) <- c("variable","Group")
    ARG_type3 <- merge(ARG_type2,group)

    ARG_type4 <- spread(ARG_type3,Type,value)
    ARG_type4 <- ARG_type4[,2:ncol(ARG_type4)]
    type.sample <- colnames(ARG_type4[2:ncol(ARG_type4)])

    ARG_type3 <- ARG_type3 %>%
        group_by(Type,Group) %>%
        summarise(Value = mean(value))
    ARG_type3 <- spread(ARG_type3,Type,Value)
    dd <- matrix(data = NA,nrow = nrow(ARG_type3),ncol = ncol(ARG_type3)-1)
    colnames(dd) <- colnames(ARG_type3)[2:ncol(ARG_type3)]
    dd <- as.data.frame(dd)
    for (i in 1:nrow(ARG_type3)) {
        dd[i,] <- ARG_type3[i,2:ncol(ARG_type3)]/ARG_type3[1,2:ncol(ARG_type3)]
    }
    ARG_type5 <- dd
    ARG_type5$Group <- ARG_type3$Group
    ARG_type5 <- ARG_type5[,c("Group",colnames(ARG_type5)[1:(ncol(ARG_type5)-1)])]
    ARG_type6 <- melt(ARG_type5)
    ARG_type6$value2 <- log2(ARG_type6$value)

    ARG_type_abun <- ggplot(ARG_type6,aes(Group,variable)) +
        geom_tile(aes(fill = value2)) +
        geom_text(aes(label = sprintf("%0.2f",round(value,digits = 2)))) +
        scale_fill_gradient2(low = "blue",mid = "grey90",high = "Red",midpoint = 0)+
        theme_bw()+
        theme(panel.grid=element_blank(),
              panel.background = element_rect(fill = "transparent",colour = NA),
              axis.ticks = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y = element_blank(),
              axis.text.y=element_blank(),
              axis.text.x=element_blank(),
              legend.position = "none",
              plot.margin = unit(c(0,0,0,0),"cm"))

    ARG_type_abun2 <- ggplot(ARG_type6,aes(Group,variable)) +
        geom_tile(aes(fill = value2)) +
        geom_text(aes(label = sprintf("%0.2f",round(value,digits = 2)))) +
        scale_fill_gradient2(low = "blue",mid = "grey90",high = "Red",
                             midpoint = 0)+
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"),
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y = element_blank(),
              axis.text.y=element_text(colour='black',size=8),
              axis.text.x=element_text(colour = "black",size = 12,face = "bold",
                                       angle = 45,hjust = 1,vjust = 1),
              legend.position = "none")
    result <- list(ARG_type5,ARG_type_abun,ARG_type_abun2)
    return(result)
}
