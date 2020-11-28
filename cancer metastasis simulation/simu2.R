---
  title: ""
author: ""
output: 
  pdf_document:
  include:
  in_header: header.tex
latex_engine: xelatex
---
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=F,warning = F,message = F)
```

# qq-plot

```{r fig.height=3}
library(ggplot2)
library(reshape2)
vafline=function(vaf){
  f=seq(0.24,0.12,-0.001)
  vaf_plot=matrix(0,nrow=length(f),ncol=length(vaf))
  for (i in 1:length(vaf)){
    vaf_plot[,i]=sapply(f, FUN = function(x) sum(vaf[[i]] > x))
    vaf_plot[,i]=vaf_plot[,i]-vaf_plot[1,i]
  }
  vaf_plot=data.frame(f,vaf_plot)
  colnames(vaf_plot)=c('f','my_vaf','my_vaf_meta','my_vaf_prim')
  vaf_plot$inv_f <- (1/vaf_plot$f - 1/0.24)
  vaf_plot=melt(vaf_plot[,-1],id.vars='inv_f')
  colnames(vaf_plot)=c('inv_f','data','mf')
  return(vaf_plot)
}
setwd('C:/Users/22367/Desktop/Hu Lab Journal Club 2018-2019/2019-7-3随机模拟修改')
my_vaf_prim=as.numeric(scan('vaf_prim.txt'))
my_vaf_meta=as.numeric(scan('vaf_meta.txt'))
my_vaff_prim=as.numeric(scan('vaff_prim.txt'))
my_vaff_meta=as.numeric(scan('vaff_meta.txt'))
my_vaf=as.numeric(scan('vaf_all.txt'))
vaf=list(my_vaf,my_vaf_meta,my_vaf_prim)
vaf_data<-data.frame(my_vaf,my_vaff_meta,my_vaff_prim)
vaf_plot=vafline(vaf)
ggplot(vaf_plot,aes_string(x='inv_f',y='mf'))+geom_point(pch=20)+
  facet_wrap(~data,nrow=1,labeller=label_both)+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(),strip.background = element_rect(fill="white"))

```

# histogram

```{r fig.height=14}
par(mfrow=c(3,2))
hist(my_vaf[my_vaf>0],breaks=seq(0,1,0.000001),xlab='allelic frequency f (0,0.0001)',main='my_vaf',col='blue',xlim=c(0,0.0001))
hist(my_vaf[my_vaf>0],breaks = 500,xlab = 'allelic frequency f (0,0.3)',main='my_vaf',col='blue',xlim=c(0,0.3))
hist(my_vaf_prim[my_vaf_prim>0],breaks=seq(0,1,0.000001),xlab='allelic frequency f (0,0.0001)',main='my_vaf_prim',col='blue',xlim=c(0,0.0001))
hist(my_vaf_prim[my_vaf_prim>0],breaks = 500,xlab = 'allelic frequency f (0,0.3)',main='my_vaf_prim',col='blue',xlim=c(0,0.3))
hist(my_vaf_meta[my_vaf_meta>0],breaks = 500,xlab = 'allelic frequency f (0,0.3)',main='my_vaf_meta',col='blue',xlim=c(0,0.3))
hist(my_vaf_meta[my_vaf_meta>0],breaks=seq(0,1,0.00001),xlab='allelic frequency f (0,0.1)',main='my_vaf_meta',col='blue',xlim=c(0,0.1))
```

# 2d distribution

```{r echo=FALSE, fig.height=8, fig.width=8}
ggplot(vaf_data,aes(my_vaff_meta,my_vaff_prim)) +
  geom_bin2d(bins = 50)+
  geom_hline(yintercept = mean(vaf_data$my_vaff_prim,na.rm = T),colour="white",size=1)+
  geom_vline(xintercept = mean(vaf_data$my_vaff_meta,na.rm = T),colour="white",size=1)+
  labs(fill="COUNT",
       title="2-D distribution of Prim and Meta",
       x="meta",
       y="prim")

```


# 晚转移
## 参数:

- 细胞总数: 100000
- 出生率: b=log(2)
- 死亡率: d=0 （假设无死亡）
- 适合度: s=0 （完全中性进化）
- 突变率: mu=20
- 变异时刻: tevent=8

## 结果:
- 生长总时间：t约等于16
- 总变异数：约400万
- 出现转移有417个细胞；最终有133个转移的细胞,里面有5651个变异

## 其他:
- 不考虑测序，使用“真实”数据
- 假设无原始突变（no clones）

图1从左到右依次为模拟的全部细胞(两个尺度）、metasis细胞、primary细胞的M(f)~1/f的散点图；图2从上到下依次为模拟的全部细胞、metasis细胞、primary细胞的vaf直方图。
                 
                 
                 
                 
                 
                 
                 
                 