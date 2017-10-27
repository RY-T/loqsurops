do_something <- function(data_path, out_path) {
  d=read.csv(file =data_path,sep = '\t',header = FALSE)
  d$V3<-d$V1
  d$V3<-substring(d$V3,0,4)
  d$V4<-d$V1
  d$V4<-substring(d$V4,5)
  d$V4<-gsub("[^0-9]", "",d$V4)
  d <- d[-c(1)]
  
  attach(d)
  d$V5[V4==9]<-'hpRNA'
  d$V5[V4>=10 & V4<=11]<-'retrotransposon'
  d$V5[V4==12]<-'DNA Transposon'
  d$V5[is.element(V4,13:19)] <- 'Others'
  detach(d)
  smaller <- subset(d,d$V5!='retrotransposon')
  smaller<-subset(smaller,smaller$V3!='B145')
  smaller
  collapse=aggregate(V2~V5+V3,smaller,sum)
  png(filename = out_path)
  plot<-barplot(collapse$V2, names.arg = collapse$V3, col=c(2:4), legend=unique(collapse$V5))
  dev.off()
}

do_something(snakemake@input[[1]], snakemake@output[[1]])