# required libraries
library("ggplot2")
library("latex2exp")

# Define custom locations
means <- c(-2.1,-1.9,-0.1, 0.1, 1.9, 2.2)

# Define functions to plot
true_dens <- function(x) {
  0.3*dnorm(x,-2,0.5) + 0.5*dnorm(x, 0,0.5) + 0.2*dnorm(x, 2,0.5)
}
dens_w1 <- function(x) {
  0.15*dnorm(x,means[1],0.5) + 0.15*dnorm(x,means[2],0.5) +
  0.25*dnorm(x,means[3],0.5) + 0.25*dnorm(x,means[4],0.5) +
  0.1*dnorm(x,means[5],0.5) + 0.1*dnorm(x,means[6],0.5)
}
dens_w2 <- function(x) {
  0.29*dnorm(x,means[1],0.5) + 0.01*dnorm(x,means[2],0.5) +
  0.01*dnorm(x,means[3],0.5) + 0.49*dnorm(x,means[4],0.5) +
  0.199*dnorm(x,means[5],0.5) + 0.001*dnorm(x,means[6],0.5)
}

# Generate plot
plt_counterexample <- ggplot() + xlim(c(-4,4)) +
  geom_function(fun = true_dens, aes(color = "A"), lty=2, n = 1000) +
  geom_function(fun = dens_w1, aes(color = "B"), n = 1000) +
  geom_function(fun = dens_w2, aes(color = "C"), n = 1000) +
  scale_color_manual(values = c("A"='black', "B"='steelblue', "C"='darkorange'),
                     labels = c("A"="True", "B"=unname(TeX("$w_{\\,1}$")), "C"=unname(TeX("$w_{\\,2}$"))),
                     guide = guide_legend(direction = "horizontal", override.aes = list(linetype = c(2,1,1)))) +
  geom_point(data = data.frame(x = means, y = rep(0,6)), aes(x=x,y=y), pch=16, color='darkred') +
  # ylab(unname(TeX("$f_0(x)$"))) + xlab(unname(TeX("x"))) +
  ylab(NULL) + xlab(NULL) +
  theme(legend.title = element_blank(), legend.position = "none")

# Show
x11(height = 2, width = 4); plt_counterexample

# Save
pdf("curves_counterExample.pdf", height = 2, width = 4); plt_counterexample; dev.off()
