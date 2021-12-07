# Manuscript:   Multi-omic machine learning predictor of breast cancer therapy response
# Author:       Stephen-John Sammut
# Description:  ggplot theme for manuscript
#=========================================================================

library (ggplot2)
library (ggpubr)

theme_manuscript <- function(base_size=12, base_family="arial") {
  library(grid)
  (ggthemes::theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",hjust = 0.5,size = base_size),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title.y = element_text(angle=90,vjust =2,size = base_size),
            axis.title.x = element_text(vjust = -0.2,size = base_size),
            axis.text = element_text(size = base_size), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(size=base_size)
    ))
}

scale_fill_RCB <- function(...){
  library(scales)
  discrete_scale("fill","sjs",manual_pal(values = c("#20A39E","#ffe671","#fdb462","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_fill_pCR_RD <- function(...){
  library(scales)
  discrete_scale("fill","sjs",manual_pal(values = c("#20A39E","#fdb462")), ...)
}

scale_colour_RCB <- function(...){
  library(scales)
  discrete_scale("colour","sjs",manual_pal(values = c("#47A8BD","#ffe671","#fdb462","#ef3b2c","#662506","#5AB4E5","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_pCR_RD <- function(...){
  library(scales)
  discrete_scale("colour","sjs",manual_pal(values = c("#20A39E","#fdb462")), ...)
}


