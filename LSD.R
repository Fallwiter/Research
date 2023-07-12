library(tidyverse)
library(broom)
library(agricolae)

rm(list=ls())
#导入数据
exdat <- read.table("alpha12.txt", header = T, sep = "\t")
str(exdat)

#数据转换
Ldat <- pivot_longer(exdat, cols = -group, 
                    names_to = c("value"), 
                    values_to = "income")

#方差齐性检验（这里选择 bartlett.test）
bartst <- Ldat %>%
  group_by(value) %>% 
  nest () %>% 
  mutate(model = map(data, ~ bartlett.test(income ~ group,data = .))) %>%
  mutate(res = map(model, ~tidy(.))) %>% 
  unnest(res) %>%
  dplyr::select(p.value) 
print(bartst)#p.value >0.05表示方差齐整可以进行单因素方差分析

#批量单因素方差分析
slist<- unique(Ldat$value ) %>% as.vector()
result<- data.frame()

for(i in 1:length(slist)){
  dfa <- Ldat %>%
    filter(value == slist[i])#拆分数据集
  
  tuk <- aov(income ~ group, dfa)
  lsd <- LSD.test(tuk, "group", p.adj = "BH") #使用LSD.test
  
  lsdr <- lsd$groups %>% as.data.frame() %>%
    mutate(group= rownames(lsd$groups),
           value= slist[i]) %>%
    dplyr::select(-income) 
  
  result<- rbind(result, lsdr)
}

#计算平均值和标准差
dfsa <- Ldat %>% group_by(group,value) %>%
  summarize(mean = mean(income),
            sd = sd(income),              
            #se = sd(income) / sqrt(length(income))
            )

#合并数据框
dfsp <- merge(dfsa, result, by=c("group", "value")) %>%
  arrange(by=value)


