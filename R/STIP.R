#' Perform STIP
#'
#' State Transition Inference Prediction
#'
#' This function generates a plot that performs (STIP) State Transition Inference Prediction.
#' @param fit a numeric matrix of gene expression values. Each row is a gene and each column is a cell. Cells are ordered according to pseudotime.
#' @param gl a character vector of gene names that will be highlighted.
#' @return A ggplot object.
#' @export
#' @import reshape2 gplots ggplot2 ggExtra
#' @author Zhicheng Ji, Zeyu Chen <zji4@@zji4.edu>
#' @examples
#' set.seed(12345)
#' fit <- t(apply(matrix(rnorm(1000),nrow=10),1,sort))
#' row.names(fit) <- LETTERS[1:10]
#' gl <- c("C", "G", "E", "I", "A")
#' STIP(fit,gl)

STIP <- function(fit,gl) {
      dn <- dimnames(fit)
      fit <- t(apply(fit,1,scale))
      dimnames(fit) <- dn
      gene <- row.names(fit)
      
      zpdirection <- fit[,1] < fit[,ncol(fit)]
      
      zp <- apply(fit,1,function(sf) {
            which(sapply(1:(length(sf)-1),function(i) sf[i]*sf[i+1] < 0))
      })
      zpnum <- sapply(zp,length)
      inczp <- names(which(zpdirection[zpnum==1]))
      deczp <- names(which(!zpdirection[zpnum==1]))
      multipoint <- names(zpnum)[zpnum > 1]
      
      geneorder <- NULL
      if (length(deczp) > 0) {
            geneorder <- c(geneorder,names(sort(unlist(zp[deczp]),decreasing = F)))
      }
      if (length(multipoint) > 0) {
            geneorder <- c(geneorder,names(sort(sapply(zp[multipoint],function(i) i[1]))))
      }
      if (length(inczp) > 0) {
            geneorder <- c(geneorder,names(sort(unlist(zp[inczp]))))
      }
      
      geneorder <- rev(geneorder)
      plotdata <- fit[geneorder,]
      plotdata <- melt(plotdata)
      colnames(plotdata) <- c("Gene","Pseudotime","Expression")
      
      p1 <- ggplot(plotdata,aes(Pseudotime,Gene,fill=Expression)) + geom_tile() + theme_classic() + scale_fill_gradient2(low="blue",high="red",midpoint=0) + theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.line.y = element_blank(),plot.margin=margin(2,0,2,0)) + scale_x_continuous(expand=c(0.001,0.001))
      
      yax <- rep("A",length(gene))
      inid <- which(gl %in% gene)
      oid <- which(!gl %in% gene)
      yaxglid <- round(seq(1,length(gene)-2,length.out=length(inid)))
      yax[yaxglid] <- gl
      yax[setdiff(1:length(yax),yaxglid)] <- setdiff(gene,gl)
      
      p2 <- ggplot() + geom_point(data=data.frame(gene=factor(yax,levels=yax),x=1),aes(x=x,y=gene),col="white") + geom_text(data=data.frame(text=gl,id=gl,x=1),aes(x=x,y=id,label=gl),size=3) + geom_segment(data=data.frame(x=1.5,xend=2,y=gl,yend=yax[match(gl,geneorder)]),aes(x=x,y=y,xend=xend,yend=yend),size=0.1) + theme_classic() + coord_cartesian(xlim=c(0.5,2)) + theme(axis.title = element_text(color="white"),axis.line = element_line(color="white"),axis.text = element_text(color="white"),axis.ticks = element_line(color="white"),plot.margin=margin(2,0,2,2)) + scale_x_continuous(expand=c(0,0))
      
      gridExtra::grid.arrange(p2,p1,nrow=1,layout_matrix=cbind(1,2,2,2))
}
