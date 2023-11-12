#' Draw a heatmap and allow user to customize the title
#' @param host Host abundance of functional gene category
#' @param title Custom title for heat map
#' @examples
##' data(ko)
##' data(group)
##' data(Gene)
##' data(tax)
##' data(abundance)
#'
#' result <- host.ratio(AnCF[[1]],title)
#'
#' @author contact: Zhengfu Yue \email{yuezhengfu2011@163.com} Liang Zeng \email{zengliang@biozeron.com}
#' @export

host.ratio <- function(host,title){
    if (nrow(host) > 6) {
        host <- host[1:5,]
    }
    pheatmap(host,fontsize = 10,cluster_rows = FALSE,fontface = "bold",cluster_cols = FALSE,
             cellwidth = 30,cellheight = 20,legend = FALSE,breaks = c(seq(0,1,by = 0.01)),
             color = colorRampPalette(c("white","Red"))(100),border_color = "black")
    grid.text(title,hjust = 0,x = 0.04,y = 0.92,gp = gpar(font = 2,size = 1.2))
}
