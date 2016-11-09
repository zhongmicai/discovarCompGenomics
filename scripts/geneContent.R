install.packages('Hmisc')
library(Hmisc)
library(ggplot2)
library(reshape2)

buscos <- read.csv("data/BUSCO2_Davey.csv")
buscosm <- melt(cbind(buscos, ind = buscos$X), id.vars = c('ind'))
buscosm <- subset(buscosm, variable!="X" & variable!="Duplicated..")
buscosm$value=as.numeric(buscosm$value)
buscosm$ind <- factor(buscosm$ind, levels = c("Agraulis vanillae","Dryas iulia",
"Eueides tales","Heliconius besckei","Heliconius burneyi","Heliconius cydno","Heliconius demeter" ,
"Heliconius elevatus","Heliconius erato mother","Heliconius erato x himera F1","Heliconius hecale",
"Heliconius hecale old","Heliconius hecalesia" ,"Heliconius himera","Heliconius himera father",
"Heliconius melpomene","Heliconius numata","Heliconius pardalinus","Heliconius sara","Heliconius telesiphe",
"Heliconius telesiphe (contaminated)","Heliconius timareta","Laparus doris","Heliconius melpomene Hmel1.1",
"Heliconius melpomene Hmel2","Heliconius erato v1","Bicyclus anynana v1.0","Bombyx mori GCA_000151625.1",
"Chilo suppressalis CsuOGS1.0","Danaus plexippus v3","Lerema accius v1.1","Manduca sexta Msex_1.0" ,
"Melitaea cinxia MelCinx1.0","Neruda aoede (contaminated)","Papilio glaucus v1.1",
"Papilio polytes Ppol_1.0","Papulio xuthus Pxut_1.0","Pieris napi DAS5","Plodia interpunctella v1",
"Plutella xylostella DBM_FJ_v1.1"))
buscosm <- subset(buscosm, ind %nin% c("Heliconius telesiphe (contaminated)","Heliconius hecale old","Neruda aoede (contaminated)"))


ggplot(buscosm) + 
  geom_bar(aes(x = ind, y = value,fill=variable),stat = "identity")+
  coord_cartesian(ylim = c(50, 100))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=15))+
  theme(plot.margin=unit(c(1.25,1.25,1.25,1.25),"cm"))
