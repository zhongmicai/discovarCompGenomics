library(ggplot2)

summaries <- read.delim("data/earlham_scaffoldsOnly.tsv")

plot <- ggplot(summaries)+
  geom_point(aes(x=X.Seqs,y=n50))+
  geom_text(aes(label=Species,x=X.Seqs,y=n50, vjust=-.5, hjust=-.05))+
  scale_x_continuous(limits = c(10000,73000))+
  scale_y_continuous(limits=c(15000,110000))+
  labs(x="Number of Scaffolds")
plot
