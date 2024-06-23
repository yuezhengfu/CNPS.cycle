#' Heat maps and statistical tests are drawn by group
#' @param ARG_type C.abundance, N.abundance, P.abundance, S.abundance
#' @param group Sample grouping file
#' @param Group_numb group number
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#'
#' result <- abun.heatmap.g(C.abundance,group,Group_numb)
#'
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export
#'

abun.heatmap.g <- function(ARG_type,group,Group_numb){
    ARG_type2 <- ARG_type[which(rowSums(ARG_type[,2:ncol(ARG_type)])>0),]
    ARG_type2 <- melt(ARG_type2)
    colnames(ARG_type2) <- c("Type","variable","value")
    colnames(group) <- c("variable","Group")
    ARG_type3 <- merge(ARG_type2,group)

    ARG_type4 <- spread(ARG_type3,Type,value)
    ARG_type4 <- ARG_type4[,2:ncol(ARG_type4)]
    type.sample <- colnames(ARG_type4[2:ncol(ARG_type4)])

    shapiro <- c()
    for (i in type.sample) {
        shapiro_test <- shapiro.test(ARG_type4[[i]])
        normality <- ifelse(shapiro_test$p.value > 0.05, "normality", "non-normality")
        shapiro <- c(shapiro,normality)
    }

    normality.p.value <- c()
    normality.method <- c()
    if (length(levels(ARG_type4$Group)) == 2) {
        for (i in type.sample) {
            fit1 <- t.test(as.formula(sprintf("`%s` ~ Group",i)),
                           data = ARG_type4)
            normality.p.value <- c(normality.p.value,fit1$p.value)
            normality.method <- c(normality.method, "T_test")
        }
    }else{
        for (i in type.sample) {
            fit1 <- aov(as.formula(sprintf("`%s` ~ Group",i)),
                        data = ARG_type4)
            normality.p.value <- c(normality.p.value,summary(fit1)[[1]][["Pr(>F)"]][[1]])
            normality.method <- c(normality.method, "ANOVA")
        }
    }

    non_normal.p.value <- c()
    non_normal.method <- c()
    if (length(levels(ARG_type4$Group)) == 2) {
        for (i in type.sample) {
            fit1 <- wilcox.test(as.formula(sprintf("`%s` ~ Group", i)), data = ARG_type4)
            non_normal.p.value <- c(non_normal.p.value,fit1$p.value)
            non_normal.method <- c(non_normal.method, "Wilcoxon")
        }
    }else{
        for (i in type.sample) {
            fit1 <- aov(as.formula(sprintf("`%s` ~ Group",i)),
                        data = ARG_type4)
            non_normal.p.value <- c(non_normal.p.value,summary(fit1)[[1]][["Pr(>F)"]][[1]])
            non_normal.method <- c(non_normal.method, "Kruskal-Wallis")
        }
    }

    p.value <- as.data.frame(cbind(type.sample,shapiro,normality.p.value,normality.method,non_normal.p.value,non_normal.method))
    colnames(p.value) <- c("Type","Normality.test","Normality.p.value","Normality.method","Non.normality.p.value","Non.normality.method")
    p.value$Normality.p.value <- as.numeric(p.value$Normality.p.value)
    p.value$Non.normality.p.value <- as.numeric(p.value$Non.normality.p.value)
    p.value$value <- ifelse(p.value$Normality.p.value < 0.05,1,NA)

    ARG_type3 <- ARG_type3 %>%
        group_by(Type,Group) %>%
        summarise(Value = mean(value))
    ARG_type3$Value <- log10(ARG_type3$Value)
    a <- floor(min(ARG_type3$Value[ARG_type3$Value != -Inf]))
    b <- floor(max(ARG_type3$Value[ARG_type3$Value != -Inf]))
    ARG_type3$Value <- ARG_type3$Value - a
    ARG_type3[ARG_type3==-Inf] <- 0
    if (b > 0) {
        i <- seq(a+1,b+1,1)
    }else{
        i <- seq(a+1,b,1)
    }

    if (b > 0) {
        breaks = c(0:(floor(max(ARG_type3$Value))+1))
    }else{
        breaks = c(0:floor(max(ARG_type3$Value)))
    }
    ARG_type_abun <- ggplot(ARG_type3,aes(Group,Type)) +
        geom_tile(aes(fill = Value),colour = "white") +
        geom_point(data = p.value,aes(x = Group_numb + 1,y = Type,size = value),
                   color = "#B2182B",show.legend = FALSE) +
        geom_point(data = p.value,aes(x = Group_numb + 1.5,y = Type),color = "white",
                   show.legend = FALSE) +
        scale_fill_gradientn(name = "Abundance",
                             colours = colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100),
                             breaks = breaks,
                             labels = c("N.D.",10^i)) +
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
              legend.title = element_text(color = "black",face = "bold",size = 12),
              legend.text = element_text(color = "black",size = 12,face = "bold")) +
        guides(fill = guide_colorbar(barheight = 10))
    result <- list(p.value,ARG_type_abun)
    return(result)
}
