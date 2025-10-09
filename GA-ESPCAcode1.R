setwd("E:\\R-language")
library(data.table)
gene_new_o <- as.matrix(fread("secondpureraw.txt", header = FALSE))
result.1_p2 <- read.table("secondbianguanxi.txt", quote = "\"", comment.char = "",header = FALSE)
result.1_p2<-as.matrix(result.1_p2)
anova_f<-read.csv("anova-f.csv",header = TRUE)
source('fun_SPCA.R')
n = 500
myString <- "start!"
print(myString)
edges = list()
gene1 = data.matrix(gene_new_o)
gene = t(gene1)
for (i in 1:12245300){
  edges[[i]] = c((result.1_p2[i,1]),(result.1_p2[i,2]))
}
weight_f <-as.numeric(as.character(anova_f[,1]))
w_f = weight_f
source('GA-ESPCA source code 2.R')
out3 =  ESPCA(gene, k = 5, edges, k.group=4000,we=0.6,t = 0.1,w_f = weight_f)
myString <- "end!"
print(myString)

setwd("E:\\R-language")
# If you want to save U matrix to csv file, use this line instead
# write.csv(out3$U, file = "u_matrix.csv", row.names = FALSE)
write.table(out3[["V"]], "V.txt",row.names = FALSE)
