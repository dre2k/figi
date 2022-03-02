# try using grid.arrange, see how it turns out 
library(png)
library(grid)
library(gridExtra)

a <- readPNG('~/tmp/figure1_a.png')
b <- readPNG('~/tmp/figure1_b.png')
c <- readPNG('~/tmp/figure1_c.png')
d <- readPNG('~/tmp/figure1_d.png')


grid.arrange(a,b)
grid.arrange(rasterGrob(a),rasterGrob(c),
             rasterGrob(b),rasterGrob(d),ncol=2, labels = c("A)", "B)", "C)", "D)"))

marrangeGrob(rasterGrob(a),rasterGrob(c),
             rasterGrob(b),rasterGrob(d),ncol=2, nrow = 2)

ggarrange(rasterGrob(a),rasterGrob(c),rasterGrob(b),rasterGrob(d), ncol = 2, nrow = 2, labels = c("A)", "B)", "C)", "D)"))


ggsave(file = "~/Dropbox/Working/test.png", width = 10, height = 7.5)



myplot1 <- arrangeGrob(rasterGrob(a), top = textGrob("A)", x = unit(0, "npc"), y   = unit(0.5, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=14, fontfamily="Times Roman")))
myplot2 <- arrangeGrob(rasterGrob(c), top = textGrob("B)", x = unit(0, "npc"), y   = unit(0.5, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=14, fontfamily="Times Roman")))
myplot3 <- arrangeGrob(rasterGrob(b), top = textGrob("C)", x = unit(0, "npc"), y   = unit(0.5, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=14, fontfamily="Times Roman")))
myplot4 <- arrangeGrob(rasterGrob(d), top = textGrob("D)", x = unit(0, "npc"), y   = unit(0.5, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=14, fontfamily="Times Roman")))

grid.arrange(myplot1, myplot2, myplot3, myplot4, ncol = 2)
pout <- arrangeGrob(myplot1, myplot2, myplot3, myplot4, ncol = 2)

ggsave(pout, file = "~/Dropbox/Working/test.png", width = 10, height = 7.5)

marrangeGrob(P1, P2, P3, labels = c("A", "B", "C"))
          ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
             