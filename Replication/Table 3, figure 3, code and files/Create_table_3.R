
library(dplyr)
library(countrycode)

setwd("./")
data <- read.csv("../../data/river_graph.csv")


p = 1.2
discharge_weight = 0.5
data = data[data$id_Country>0,]
data$score = data$concentration^p *data$flow_HR^discharge_weight * data$cell_distance
data$weight = data$cell_distance * data$flow_HR^discharge_weight
data$no_cont = data$cell_distance * as.integer(data$concentration==0)
data$much_cont = data$cell_distance * as.integer(data$concentration>2)
data$some_cont = data$cell_distance - data$much_cont - data$no_cont
data = data[,-1]
data_sum = aggregate(data[,-16], by=list(Category=data$Country), FUN=sum)
data_sum$score = data_sum$score / data_sum$weight
data_sum$score = data_sum$score/mean(data_sum$score)
data_sum$cell_distance = data_sum$cell_distance/1000
data_sum$no_cont = data_sum$no_cont/1000
data_sum$much_cont = data_sum$much_cont/1000
data_sum$some_cont = data_sum$some_cont/1000

data_sum$no_cont_perc = data_sum$no_cont/data_sum$cell_distance
data_sum$some_cont_perc = data_sum$some_cont/data_sum$cell_distance
data_sum$much_cont_perc = data_sum$much_cont/data_sum$cell_distance
data_sum = data_sum[,c(1, 20, 22, 23, 24, 25, 26, 27)]

# calculate ginis
for (y in 1:length(unique(data$id_Country))){
  data_country <- data[data$id_Country==y,]
  data_country <- data_country[order(data_country$concentration), ]
  
  data_country$weight_perc = cumsum(data_country$weight)/sum(data_country$weight)
  data_country$concentration_perc = cumsum(data_country$concentration*data_country$weight)/sum(data_country$concentration*data_country$weight)
  
  Gini = sum((data_country$weight/ sum(data_country$weight)) * data_country$concentration_perc)
  Gini = (0.5 - Gini)*2
  data_sum$gini[data_sum[,1] == data[data$id_Country==y,][1,16]] = Gini
}

hydroWASTE = read.csv("../../data/hydroWASTE rivers affected.csv")
data_sum$Category = countrycode(data_sum$Category, 'country.name', 'iso3c')
names = colnames(data_sum)
names[1] = "Country"
colnames(data_sum) = names
data_sum = merge(data_sum, hydroWASTE, left_join="Country", right_join="Country")
data_sum$wOtter_affected = data_sum[,4] + data_sum[,5]
data_sum$Ratio = data_sum[,11]/data_sum[,10]
colnames(data_sum) = c("Country", "Score", "0 (km)", ">2 (km)", "0-2 (km)", "0 (%)", "0-2 (%)", ">2 (%)", "Gini", "HydroWASTE (km)", "wOtter (km)", "Ratio")
data_sum = data_sum[, c(1, 2, 9, 3, 6, 5, 7, 4, 8, 11, 10, 12)]
write.csv(data_sum, './table2.csv')
