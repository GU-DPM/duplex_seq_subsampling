plotSpectrumTumorNormal <- function(inData) {
  library(ggplot2)
  myPlot <- ggplot(data=inData, 
                   aes(x=Measure, y=Value)) + 
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour="black"), 
          axis.text.x = element_text(angle=90,hjust=1)) + 
    labs(x="Mutation Type", y="Proportion of Mutations") + 
    geom_point(aes(colour=Type), 
               position=position_jitterdodge(jitter.width = .5, 
                                             dodge.width = .8)) + 
    geom_boxplot(aes(fill = Type), 
                 outlier.color = NA, position = position_dodge(width=.8)) + 
    scale_fill_manual(values = c(NA,NA))
  return(myPlot)
}
