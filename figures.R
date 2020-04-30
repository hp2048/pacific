library(plotly)
library(ggplot2)

##location accuracy for prediction along the genome
aln <- read.table("/media/dstore/covid19/sarspredictions.txt", header = F, sep = "\t")
colnames(aln) <- c("position","ptype","counts")
aln %>% group_by(position) %>% 
  arrange(position,ptype) %>% 
  mutate(Y = counts/sum(counts)) %>% 
  filter(ptype != "TruePositive") %>% 
  mutate(z = cumsum(Y)) %>% 
  plot_ly(type = 'scatter', x = ~position, y = ~z, color = ~ptype, mode = 'lines', fill = 'tonexty')

##pairwise alignments of genome
aln <- read.table("///media/dstore/covid19/pairwise.aln.txt", header = F, sep = "\t")
aln$qg <- sapply(strsplit(as.character(aln$V1),":"), `[`, 3)
aln$tg <- sapply(strsplit(as.character(aln$V2),":"), `[`, 3)
aln$qgroup <- sapply(strsplit(as.character(aln$V1),":"), `[`, 2)
aln$tgroup <- sapply(strsplit(as.character(aln$V2),":"), `[`, 2)
h <- aln %>% filter(V1 != V2) %>% group_by(qgroup, qg,tgroup, tg) %>% summarise(avgpid = mean(V3))
plot_ly(x=paste(h$qgroup,h$qg,sep=":"), y=paste(h$tgroup,h$tg,sep=":"), z = h$avgpid, type = "heatmap")






