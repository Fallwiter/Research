#library
library(vegan)

#清除环境变量
rm(list = ls())

#OTU数据集
otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', 
                   stringsAsFactors = FALSE, check.names = FALSE) 

##检查数据集是否抽平，理论上alpha多样性使用的测序数据不应该进行抽平
##样品otu和相等及抽平
colSums(otu)#参考数据集未抽平

otu <- t(otu)
###合并计算函数
###勿动
Alpha_diversity <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Obs <-  est[1, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')   #注意，这里是Gini-Simpson 指数
  Pielou <- Shannon / log(Obs, base)            #Pielou 均匀度
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)  #反映测序深度的指标
  result <- rbind(est, Shannon, Simpson,
                  Pielou, goods_coverage)
  #计算谱系多样性
  if (!is.null(tree)) {
    Pd <- pd(x, tree, include.root = FALSE)[1]
    Pd <- t(Pd)
    result <- rbind(result, Pd)
  }
  result <- as.data.frame(t(result))
  return(result)
}

#无tree文件，不计算谱系多样性
alpha_otu <- Alpha_diversity(otu)

#输出
alpha_otu <- cbind(sample=row.names(alpha_otu),alpha_otu)
write.table(alpha_otu,file="alpha_otu.txt",sep="\t",row.names=F,quote=F)

#有tree文件
#library(picante)
#加载进化树文件
#tree <- read.tree('otu_tree.tre')
#alpha <- Alpha_diversity(otu, tree)
