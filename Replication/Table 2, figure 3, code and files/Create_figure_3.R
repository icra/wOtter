library(ggplot2)
library(extrafont)

loadfonts(device = "win")

setwd("./")

data <- read.csv("../../data/river_graph.csv")
data = data[data$id_Country>0,]
discharge_weight = 0.5
data$weight = data$cell_distance * data$flow_HR^discharge_weight

data <- data[order(data$concentration), ]
elements = length(data[,1])

data$weight_perc = cumsum(data$weight)/sum(data$weight)
data$concentration_perc = cumsum(data$concentration*data$weight)/sum(data$concentration*data$weight)

Gini = sum((data$weight/ sum(data$weight)) * data$concentration_perc)
Gini = (0.5 - Gini)*2

Gini_label = paste('G=', substr(as.character(Gini), 1, 4), sep='')


ggplot(data, aes(x=weight_perc, y=concentration_perc)) + geom_line() +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0.015)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0),position = "right") +
  labs(x="Percentage of (weighted) rivers",y="Percentage of contamination") +
  geom_abline() +
  geom_polygon(data = data, aes(x=weight_perc, y=concentration_perc),
               fill="cyan4", alpha=0.25) +
  annotate("text", x = 0.53, y=0.40, label = Gini_label, size=16,
           family="serif") + 
  theme(axis.text.x=element_text(size=25, family="serif"),
        axis.text.y=element_text(size=25, family="serif"),
        axis.title.x =element_text(size=30, family="serif"),
        axis.title.y =element_text(size=30, family="serif"))
ggsave("./lorenz_curve.png", width = 20, height = 10)