
library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(metR)

#########################
### SATF Illustration ###
#########################

geom_double_arrow <- function(p, x1, y1, x2, y2) { 
  p + geom_segment(aes(x=x1, xend=x2, y=y1, yend=y2), arrow = arrow(type = "closed", length = unit(.25,"cm")), size=.1) +
      geom_segment(aes(x=x2, xend=x1, y=y2, yend=y1), arrow = arrow(type = "closed", length = unit(.25,"cm")), size=.1)
}

dprime <- function(t, lambda=.5, beta=1.5, delta=1.5) ifelse( t >= delta, lambda * (1-exp(-(t-delta)/beta)), 0)

y1 = 0.05
y2 = 0.1

x <- seq(0, 10, 1)
p <- ggplot(data.frame(x=x, y= dprime(x) ), aes(x, y)) + stat_function(fun = dprime)

p <- p %>% geom_double_arrow(x1=0, x2=1.5, y1=y1, y2=y1) + 
            geom_label(aes(x=0.75, y=y1, label = "δ"), colour = "black",  position= position_dodge(width=1))

p <- p %>% geom_double_arrow(x1=1.5, x2=1.5+1.5, y1=y1, y2=y1) +
            geom_label(aes(x=2.25, y=y1, label = "β^{-1}"), colour = "black",  position= position_dodge(width=1), parse=T)

p <- p + geom_segment(aes(x=1.5, xend=1.5, y=0, yend=y1), color="black", linetype = "dotted", size=.1) + 
          geom_segment(aes(x=3, xend=3, y=.5*(1-exp(-1)), yend=0), color="black", linetype = "dotted", size=.1) + 
          theme_bw() +
          scale_x_continuous(breaks = seq(0, 10, 1.5))

p <- p + geom_segment(aes(x=0, xend=0, y=0, yend=y1), color="black", linetype = "dotted", size=.1)

p <- p + geom_segment(aes(x=-Inf, xend=10, y=.5, yend=.5), color="black", linetype = "dotted", size=.1) +
          geom_segment(aes(x=-Inf, xend=3, y=.5*(1-exp(-1)), yend=.5*(1-exp(-1)) ), color="black", linetype = "dotted", size=.1) +
          geom_segment(aes(x=-Inf, xend=Inf, y=0, yend=0 ), color="black", linetype = "dotted", size=.1)

p <- p + geom_label(aes(x=1, y=.5*(1-exp(-1)), label = "0.63·λ"), colour = "black",  position= position_dodge(width=1)) +
          geom_label(aes(x=1, y=.5, label = "λ"), colour = "black",  position= position_dodge(width=1))

p <- p + scale_y_continuous(breaks=c(0, .5), labels = c(0,"λ") ) + scale_x_continuous(breaks=NULL)
p <- p + theme(axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("time") + ylab("accuracy (d')")

p <- p + geom_point()

z = .8
ggsave(p, file = "../figures/illustrations/illustration_satf.pdf", height = z*3.5, width = z*7.26, device = cairo_pdf)

