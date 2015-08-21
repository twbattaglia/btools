plot_nice_boxplot <- function(phylobj, x, color, measures = "Observed", save = T, name = "Figure.pdf"){
  
  # Make plot
  plot_data <- plot_richness(physeq = phylobj, 
                             x = x, 
                             measures = measures, 
                             color = color)
  p <- plot_data + geom_boxplot(aes(fill = Treatment), outlier.colour = NA, alpha = 0.80, position = position_dodge(width = .9)) + theme_minimal(base_size = 20) +  theme(legend.position = "bottom", panel.grid.major.x = element_line(colour = "grey", size = 0.55)) + xlab("") + ylab("") + guides(fill = guide_legend(title = "Treatment")) + scale_fill_manual(values = c("Control" = "black","PAT" = "#57068c")) + scale_color_manual(values = c("Control" = "black","PAT" = "#57068c")) + geom_point(size = 3.5)
  
  if (save == TRUE){
    p <- p + ggsave(filename = name, width = 11.69, height = 8.27)
    return(p)
  }
  else{
    return(p)
  }
}



plot_nice_boxplot_table <- function(phylobj, x, group, measures = "Observed", test_type = "nonparametric", save = T, name = "Figure.pdf"){
  
  # Library 
  library(gridExtra)
  
  # Get Pvalue table
  pvalue_table <- compare_alpha_diversity(physeq = phylobj, 
                                          x = x, 
                                          group = group, 
                                          diversity = measures, 
                                          test_type = test_type)
  
  # Make boxplot of data
  boxplot <- plot_nice_boxplot(phylobj, 
                               x = x, 
                               measures = measures, 
                               color = group, 
                               save = F)
  
  # Combine together
  pdf(file = name, paper = "a4r", width = 11.69, height = 8.27)
  grid.arrange(tableGrob(d = pvalue_table, 
                         show.rownames = F, 
                         theme = theme.white()), 
               boxplot, 
               nrow = 2, 
               as.table = TRUE)
  dev.off()
  
  return(grid.arrange(tableGrob(d = pvalue_table, 
                                show.rownames = F, 
                                theme = theme.white()), 
                      boxplot, 
                      nrow = 2, 
                      as.table = TRUE))
  
}