library(ggplot2)
library(reshape2)

buscos <- read.csv("data/BUSCO2_Davey.csv")
buscosm <- melt(cbind(buscos, ind = buscos$X), id.vars = c('ind'))
buscosm <- subset(buscosm, variable!="X" & variable!="Duplicated..")
buscosm$value=as.numeric(buscosm$value)

ggplot(buscosm) + 
  geom_bar(aes(x = ind, y = value,fill=variable),stat = "identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
