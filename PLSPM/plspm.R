##plspm需要通过github安装
#library(devtools)
#install_github("gastonstat/plspm")
library(plspm)
library(dplyr)

#清除环境
rm(list=ls())

exdat <- read.table("data.txt", header = T, sep = "\t")

#查看数据类型
str(exdat)

#将“chr”更改为“factor”
exdat$sample <- factor(exdat$sample)
exdat$group <- factor(exdat$group)
str(exdat)

#通过 0-1 矩阵描述变量之间的关联，其中 0 代表变量间没有关联，1 代表有关联
#构建模型矩阵，上半矩阵填充为0，有关联为1(不清楚可以去查询plspm包的帮助文档)
Biocarbon <- c(0, 0, 0, 0, 0, 0, 0)
Water <- c(0, 0, 0, 0, 0, 0, 0)
N20 <- c(1, 1, 0, 0, 0, 0, 0)
Soil <-c(1, 1, 1, 0, 0, 0, 0)
MBN <- c(1, 1, 1, 0, 0, 0, 0)
Nitrifiers <- c(1, 1, 1, 1, 1, 0, 0)
Denitrifiers <- c(1, 1, 1, 1, 1, 0, 0)

#矩阵(行名等于列名，并且按照上述顺序绑定)
edu_path <- rbind(Biocarbon, Water, N20, Soil, MBN,Nitrifiers, Denitrifiers)
colnames(edu_path) = rownames(edu_path)

#绘制内部矩阵
innerplot(edu_path, box.size = 0.1)

#构建外部模型（指定潜在变量），根据研究的相同类型对进行绑定
##举例说明为Biocarbon对应数据框中的第一列，Soil对应数据框中的第四到七列
edu_blocks <- list(1,2,3,c(4:7),8,c(9:10),c(11:13))

#按照顺序列出每个变量的类型
#变量类型分为"A","B"两种，具体可查询plspm包的帮助文档
edu_modes <- rep("A",7)
#edu_modes <- c("A","A","A","A","A","A","B")

#plspm模型初构，boot.val=T执行引导验证，br重采样次数
edu_val <- plspm(exdat, edu_path, edu_blocks, modes = edu_modes, 
                 boot.val = TRUE, br = 200)

##检查模型一维性，需严格满足下列条件，不符合模型需重构
edu_val$unidim
#C.alpha, DG.rho >0.7
#eig.1st >1
#eig.2nd <1
##PS: 这里的模型是不合格的未进行重构,可以根据模型检测结果重构
plot(edu_val, what = "loadings")

#检测外部模型（潜在变量与变量之间的关联性）
edu_val$outer_model
#loadings >0.7, communality >0.5

#每个指标必须在自己的因子下的crossloading值最大
edu_val$crossloadings

#可视化
library(ggplot2)
ggplot(data = edu_val$outer_model,
       aes(x = name, y = loading, fill = block)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  # threshold line (to peek acceptable loadings above 0.7)
  geom_hline(yintercept = 0.7, color = 'gray50') +
  # add title
  ggtitle("Barchart of Loadings") +
  # rotate x-axis names
  theme(axis.text.x = element_text(angle = 90))


#模型整体评价 gof 指数（gof >0.7非常好）
edu_val$gof
#[1] 0.6476621
#假模型的gof也挺高的，所以说gof只能评估模型整体情况不能具体说明模型是否可靠

#模型变量的R2（Mean.Boot为重采样结果）
edu_val$boot$rsq

#路径系数检验（Mean.Boot为重采样结果）
edu_val$boot$paths

#输出模型（这个图太丑了，可以使用draw.io/）
plot(edu_val)

#变量，潜在变量之间的影响（plspm会构建变量间的影响）
path_effs <- edu_val$effects
write.table(path_effs,file="path_effs.txt",sep="\t",row.names=F,quote=F)

#输出模型所有结果
sink('plspm-summary.txt')
summary(edu_val)
sink(NULL)
