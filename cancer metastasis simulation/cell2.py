# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 21:25:35 2019
w
@author: admin
"""
import numpy as np
import math
import random
def initiate(mu,clones=0):       # 细胞生长起始函数
    mutid=1
    cell=[]
    nmut=clones
    cell.append(list(range(1,clones+1)))
    mutid=mutid+nmut
    return cell,mutid
def mutant(celli,mu,mutid):     # 细胞分裂突变函数
    nmut=int(np.random.poisson(mu,1))
    celli=celli+list(range(mutid,mutid+nmut))
    mutid=mutid+nmut
    return celli,mutid
def copycell(celli):           # 细胞分裂，母细胞一分为二
    cellj=celli
    return(cellj)
def grow(b,d,nend,mu,s,tevent,t=0):    # 细胞生长函数，可调参数：出生率、死亡率、突变率、适合度、出现亚克隆的时间、起始时刻
    br=[b]
    dr=[d]
    s=[0]+s
    numclone=len(s)-1
    if numclone>0 :     # 判断是否加入亚克隆
        for i in range(0,numclone):
            dr.append(dr[0]*np.random.uniform())
            br.append((1 + s[i+1]) * (br[0] - dr[0]) + dr[i+1])
    subclone=0
    eventcellnum=[]
    rmax=b+d
    cell,mutid=initiate(mu=mu)    # 开始生长
    celltype=[0]    # 记录细胞类型：neutral为0，其余亚克隆按出现次序依次加1
    while len(cell)<nend:
        r=np.random.uniform(0,rmax)    # 获得r值来判断生长或死亡
        nt=len(cell)
        n=random.randint(a=0,b=nt-1)    # 随机选择一个母细胞
        if r<br[celltype[n]] :    # 细胞分裂过程
            cell.append(copycell(celli=cell[n]))
            celltype.append(list(np.where([bi==br[celltype[n]] for bi in br])[0])[0])
            cell[n],mutid=mutant(celli=cell[n],mu=mu,mutid=mutid)
            cell[-1],mutid=mutant(celli=cell[-1],mu=mu,mutid=mutid)
            t=t+1/(rmax*nt)*int(np.random.exponential(1,1))
            if subclone!=numclone :    # 判断是否还需要产生亚克隆
                if t>tevent[subclone]:    # 判断是否在该步引入亚克隆
                    subclone = subclone + 1
                    celltype[n]=subclone
                    rmax=max(br[0:subclone+1])+max(dr[0:subclone+1])
                    eventcellnum.append(nt+1)
        elif r>br[celltype[n]]+dr[celltype[n]] :    # 细胞既不分裂也不死亡
            t=t+1/(rmax*nt)*int(np.random.exponential(1,1))
        elif r>=br[celltype[n]] and r<=br[celltype[n]]+dr[celltype[n]] :    # 细胞死亡
            cell.pop(n)
            celltype.pop(n)
            t=t+1/(rmax*nt)*int(np.random.exponential(1,1))
            if len(cell)==0 :    # 判断细胞数是否为0，若是，重新生长
                t=0
                cell,mutid=initiate(mu=mu)
                celltype=[0]
    return cell,t,eventcellnum,celltype    # 返回所有细胞的突变位点、生长总时间、发生亚克隆时的细胞总数、所有细胞的细胞类型
def af(cell):    # 获得突变位点的突变频率
    mutants=sum(cell,[])
    AF=[]
    for i in range(min(mutants),max(mutants)+1):
        AF.append(mutants.count(i))
    return AF
def cal_vaf(AF,cellnum,ploidy=2,detectlimit=0.01,read_depth=100):    # 模拟测序过程产生vaf
    AF=[x/ploidy for x in AF]
    AF=list(filter(lambda x: x>detectlimit*cellnum,AF))
    depth = list(np.random.poisson(read_depth,len(AF)))
    samp_alleles = list(map(lambda x,y: np.random.binomial(x,y/cellnum),depth,AF))
    vaf=list(map(lambda x,y: x/y,samp_alleles,depth))
    return vaf     
b=math.log(2) # 出生率
d=0 # 死亡率
nend=10000 # 细胞总数
s=[0.2] # 适合度
mu=20 # 突变率
tevent=[1.5] # 产生亚克隆的时刻
read_depth=200 # 测序阅读深度
detectlimit=5/read_depth # 能检测到的位点最低频率
cell,t,ecn,ctp=grow(b=b,d=d,nend=nend,mu=mu,s=s,tevent=tevent)
AF=af(cell)
vaf=cal_vaf(AF=AF,cellnum=nend,detectlimit=detectlimit,read_depth=read_depth)
with open('vaf_subclone.txt', 'w') as fl:
    for item in vaf:
        fl.write("%s\n" % item)

## 后面是对细胞进行分类，可以无视，需要的时候再用上
cell_sub_idx=list(np.where([ctpi==1 for ctpi in ctp])[0])
cell_sub=list(cell[i] for i in cell_sub_idx)
cell_neu_idx=list(np.where([ctpi==0 for ctpi in ctp])[0])
cell_neu=list(cell[i] for i in cell_neu_idx)
AF_sub=af(cell_sub)
vaf_sub=cal_vaf(AF=AF_sub,cellnum=nend,detectlimit=detectlimit,read_depth=read_depth)
with open('vaf_subset_subclone.txt', 'w') as fl:
    for item in vaf_sub:
        fl.write("%s\n" % item)
AF_neu=af(cell_neu)
vaf_neu=cal_vaf(AF=AF_neu,cellnum=nend,detectlimit=detectlimit,read_depth=read_depth)
with open('vaf_subset_neutral.txt', 'w') as fl:
    for item in vaf_neu:
        fl.write("%s\n" % item)    
    

    
