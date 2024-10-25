plotGseaTable <- function(pathways, stats, fgseaRes,
                          gseaParam=1,
                          colwidths=c(5, 3, 0.8, 1.2, 1.2),
                          render=TRUE) {

    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    pathways <- lapply(pathways, function(p) {
        unname(as.vector(na.omit(match(p, names(statsAdj)))))
    })

    # fixes #40
    pathways <- pathways[sapply(pathways, length) > 0]

    ps <- lapply(names(pathways), function(pn) {
        p <- pathways[[pn]]
        # annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
        annotation <- fgseaRes[match(pn, fgseaRes$ID), ]  #改动和clusterprofiler列名一致
        list(
            textGrob(pn, just="right", x=unit(0.95, "npc")),
            ggplot() +
                geom_segment(aes(x=p, xend=p,
                                 y=0, yend=statsAdj[p]),
                             size=0.2) +
                scale_x_continuous(limits=c(0, length(statsAdj)),
                                   expand=c(0, 0)) +
                scale_y_continuous(limits=c(-1, 1),
                                   expand=c(0, 0)) +
                xlab(NULL) + ylab(NULL) +
                theme(panel.background = element_blank(),
                      axis.line=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank(),
                      panel.grid = element_blank(),
                      axis.title=element_blank(),
                      plot.margin = rep(unit(0,"null"),4),
                      panel.spacing = rep(unit(0,"null"),4)
                ),
            textGrob(sprintf("%.2f", annotation$NES)),
            textGrob(sprintf("%.1e", annotation$pvalue)),  #改动和clusterprofiler列名一致
            textGrob(sprintf("%.1e", annotation$p.adjust)) #改动和clusterprofiler列名一致
            )
    })


    rankPlot <-
        ggplot() +
        geom_blank() +
        scale_x_continuous(limits=c(0, length(statsAdj)),
                           expand=c(0, 0)) +
        scale_y_continuous(limits=c(-1, 1),
                           expand=c(0, 0)) +
        xlab(NULL) + ylab(NULL) +
        theme(panel.background = element_blank(),
              axis.line=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid = element_blank(),
              axis.title=element_blank(),
              plot.margin = unit(c(0,0,0.5,0), "npc"),
              panel.spacing = unit(c(0,0,0,0), "npc")
        )

    grobs <- c(
        lapply(c("Pathway", "Gene ranks", "NES", "pval", "padj"), textGrob),
        unlist(ps, recursive = FALSE),
        list(nullGrob(),
             rankPlot,
             nullGrob(),
             nullGrob(),
             nullGrob()))

    # not drawing column if corresponding colwidth is set to zero
    grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))


    p <- arrangeGrob(grobs=grobs[grobsToDraw],
                 ncol=sum(colwidths != 0),
                 widths=colwidths[colwidths != 0])

    if (render) {
        grid.draw(p)
    } else {
        p
    }
#	return(annotation)
}