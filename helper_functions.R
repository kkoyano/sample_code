#!Rscript

# Name: Kiku Koyano
# Date: 05/06/2021
# Description: common functions used to plot figures 

libraries = c("data.table", "dplyr", "optparse", "ggpubr", "tidyr", "ggplot2",
              "ggthemes", "purrr", "forcats", "RColorBrewer", "grid", "gridExtra", "readxl", "rstatix", "ggpubr") 
shh = lapply(libraries, require, quietly = T, character.only = TRUE) 

grid.ftable <- function(d, padding = unit(4, "mm"), ...) {
  
  nc <- ncol(d)
  nr <- nrow(d)
  
  ## character table with added row and column names
  extended_matrix <- cbind(c("", rownames(d)),
                           rbind(colnames(d),
                                 as.matrix(d)))
  
  ## string width and height
  w <- apply(extended_matrix, 2, strwidth, "inch")
  h <- apply(extended_matrix, 2, strheight, "inch")
  
  widths <- apply(w, 2, max)
  heights <- apply(h, 1, max)
  
  padding <- convertUnit(padding, unitTo = "in", valueOnly = TRUE)
  
  x <- cumsum(widths + padding) - 0.5 * padding
  y <- cumsum(heights + padding) - padding
  
  rg <- rectGrob(x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 width = unit(widths + padding, "in"),
                 height = unit(heights + padding, "in"))
  
  tg <- textGrob(c(t(extended_matrix)), x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 just = "center")
  
  g <- gTree(children = gList(rg, tg), ...,
             x = x, y = y, widths = widths, heights = heights)
  
  grid.draw(g)
  invisible(g)
}

ttheme_default = ttheme_default(colhead=list(bg_params=list(fill="white")),
                                core = list(bg_params = list(fill= "white")) )
             
                                                                
grid.ftable(grob, gp = gpar(fill = NA))
output_table = function(output_dir, filename, grob, width = 7, height = 8){
  # note: grobs can be either a data.frame to turn into a grob object. or it is a list grobs 
  output_fn = paste(output_dir, filename, sep = '')
  message(output_fn)
  pdf(output_fn, width = width, height = height )
  
  if (is.data.frame(grob)){
    table_grob = tableGrob(grob)
    # table_grob = tableGrob(grob, theme = ttheme_default)
    print(grid.arrange(table_grob) ) 
    
    
    
  } else if (is.list(grob)) {
    print(marrangeGrob(grob, nrow=1, ncol=1))
  } else {
    message("grob is not a dataframe or a list of grobs. please input correct value")
  }
  dev.off()
  
}

output_plots = function(output_dir, filename, plot_list, width = 5, height = 4){
  output_fn = paste(output_dir, filename, sep = '')
  pdf(output_fn, width = width, height = height )
  message(output_fn)
  print(plot_list)
  dev.off()
}
