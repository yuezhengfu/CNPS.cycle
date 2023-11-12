#' PCA analysis and visualization
#' @param ARG_sub C.abundance[,2:ncol(C.abundance)], N.abundance[,2:ncol(N.abundance)], P.abundance[,2:ncol(P.abundance)], S.abundance[,2:ncol(S.abundance)]
#' @param group Sample grouping file
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#'
#' result <- pca.arg(C.abundance[,2:ncol(C.abundance)],group)
#'
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export

pca.arg <- function(ARG_sub,group){
    colnames(group) <- c("sample","Group")
    ARG_beta <- t(ARG_sub)
    pcoa<- PCA(ARG_beta,scale.unit = FALSE,graph = FALSE)
    PC1 = pcoa$ind$coord[,1]
    PC2 = pcoa$ind$coord[,2]
    plotdata <- data.frame(rownames(pcoa$ind$coord),PC1,PC2)
    colnames(plotdata) <-c("sample","PC1","PC2")
    plotdata <- merge(plotdata,group)
    pc1 <-floor(pcoa$eig[1,2]*100)/100
    pc2 <-floor(pcoa$eig[2,2]*100)/100

    ARG_pcoa1<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title="PCA - Function genes of carbon cycling") +
        xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) +
        ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
        theme(text=element_text(size=18))+
        theme(panel.background = element_rect(fill='white', colour='black'),
              panel.grid=element_blank(),
              axis.title = element_text(color='black',size=18),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_text(colour='black', size=18),
              axis.title.y=element_text(colour='black', size=18),
              axis.text=element_text(colour='black',size=16),
              legend.title=element_text(size = 14,face = "bold"),
              legend.text=element_text(size=12),
              legend.key=element_blank(),legend.position = "right",
              legend.background = element_rect(colour = "black"))+
        theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))

    ARG_pcoa2<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        stat_ellipse(aes(fill = Group),geom = "polygon",level = 0.95,alpha = 0.3)+
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title="PCA - Function genes of cycling") +
        xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) +
        ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
        theme(text=element_text(size=18))+
        theme(panel.background = element_rect(fill='white', colour='black'),
              panel.grid=element_blank(),
              axis.title = element_text(color='black',size=18),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_text(colour='black', size=18),
              axis.title.y=element_text(colour='black', size=18),
              axis.text=element_text(colour='black',size=16),
              legend.title=element_text(size = 14,face = "bold"),
              legend.text=element_text(size=12),
              legend.key=element_blank(),legend.position = "right",
              legend.background = element_rect(colour = "black"))+
        theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))

    ARG_pcoa3<-ggplot(plotdata, aes(PC1, PC2)) +
        geom_point(aes(fill=Group),size=4.5,color = "black",shape = 21,alpha = 0.8)+
        geom_label_repel(aes(PC1,PC2,label = sample),fill = "white",color = "black",
                         box.padding = unit(0.3,"lines"),segment.colour = "grey50",
                         label.padding = unit(0.15,"lines"),size = 2) +
        geom_vline(aes(xintercept = 0),linetype="dotted")+
        geom_hline(aes(yintercept = 0),linetype="dotted")+
        scale_fill_manual(values=cbbPalette)+
        labs(title="PCA - Function genes of cycling") +
        xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) +
        ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
        theme(text=element_text(size=18))+
        theme(panel.background = element_rect(fill='white', colour='black'),
              panel.grid=element_blank(),
              axis.title = element_text(color='black',size=18),
              axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_text(colour='black', size=18),
              axis.title.y=element_text(colour='black', size=18),
              axis.text=element_text(colour='black',size=16),
              legend.title=element_text(size = 14,face = "bold"),
              legend.text=element_text(size=12),
              legend.key=element_blank(),legend.position = "right",
              legend.background = element_rect(colour = "black"))+
        theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))


    result <- list(plotdata,ARG_pcoa1,ARG_pcoa2,ARG_pcoa3)
    return(result)
}
