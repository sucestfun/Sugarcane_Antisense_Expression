# Creating a chart to visualize expression patterns among the sugarcane sense/antisense transcripts in the seven investigated experiments

# Required package
library(UpSetR)

#Getting data -------------------------------------------------------------------------------------------
data<-read.csv("sense-antisense_expression.csv", sep = ";", header = T)
data

# Making a list with all transcripts grouped by sense/antisense expression pattern
y = list(
  Concordant = data[which(data[,3]=="Concordant"),1],
  Discordant = data[which(data[,3]=="Discordant"),1],
  Neutral = data[which(data[,3]=="Neutral"),1])
y

# Plotting the data
upset(fromList(y), nsets=3, order.by = "freq", 
queries = list(list(query = intersects,
params = list ("Concordant"), color = "darkblue", active = F),
list(query = intersects,
params = list ("Discordant"), color = "darkred", active = F),
list(query = intersects,
params = list ("Neutral"), color = "khaki1", active = F),
list(query = intersects,
params = list ("Neutral", "Concordant"), color = "lightblue2", active = F),
list(query = intersects,
params = list ("Concordant", "Discordant"), color = "green4", active = F),  
list(query = intersects,
params = list ("Discordant", "Neutral"), color = "pink1")),
main.bar.color = "gray40", mainbar.y.label = "Number of SAS", 
sets.bar.color = c("darkblue", "khaki1", "darkred"), 
sets.x.label = "Set Size", point.size = 4.8, line.size = 2.0,
mb.ratio = c(0.70,0.30), 
decreasing = T, show.numbers = "yes", shade.color = "white", 
matrix.dot.alpha =0.8, attribute.plots = NULL, 
scale.intersections = "identity",
scale.sets = "identity", text.scale = 2.0, set_size.show = F
)
