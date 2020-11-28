# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 08:40:20 2019

@author: 22367
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
def grow(b,d,nend,mu,s,tevent,t=0):    # 细胞生长函数，可调参数：出生率、死亡率、突变率、适合度、出现亚克
    s=[0]+s
    numclone=len(s)-1
    br=[b]
    dr=[d]   
    if numclone>0 :     # 判断是否加入亚克隆
        for i in range(0,numclone):
            dr.append(dr[0]*np.random.uniform())
            br.append((1 + s[i+1]) * (br[0] - dr[0]) + dr[i+1])
    subclone=0
    eventcellnum=[]
    rmax=b+d
    cell,mutid=initiate(mu=mu)    # 开始生长
    celltype=[0]    # 记录细胞类型：primary为0，其余亚克隆按出现次序依次加1
    while len(cell)<nend:
        r=np.random.uniform(0,rmax)    # 获得r值来判断生长或死亡
        nt=len(cell)
        n=random.randint(a=0,b=nt-1)    # 随机选择一个母细胞
        if r<br[celltype[n]] :    # 细胞分裂过程
            cell.append(copycell(celli=cell[n]))
            celltype.append(celltype[n])
            cell[n],mutid=mutant(celli=cell[n],mu=mu,mutid=mutid)
            cell[-1],mutid=mutant(celli=cell[-1],mu=mu,mutid=mutid)
            t=t+1/(rmax*nt)*np.random.exponential(1,1)
            if subclone!=numclone :    # 判断是否还需要产生亚克隆
                if t>tevent[subclone]:    # 判断是否在该步引入亚克隆
                    subclone = subclone + 1
                    celltype[n]=subclone #让母细胞成为新亚克隆的鼻祖
                    rmax=max(br[0:subclone+1])+max(dr[0:subclone+1])
                    eventcellnum.append(nt+1)
        elif r>br[celltype[n]]+dr[celltype[n]] :    # 细胞既不分裂也不死亡
            t=t+1/(rmax*nt)*np.random.exponential(1,1)
        elif r>=br[celltype[n]] and r<=br[celltype[n]]+dr[celltype[n]] :    # 细胞死亡
            cell.pop(n)
            celltype.pop(n)
            t=t+1/(rmax*nt)*np.random.exponential(1,1)
            if len(cell)==0 :    # 判断细胞数是否为0，若是，重新生长
                t=0
                cell,mutid=initiate(mu=mu)
                celltype=[0]
    return cell,t,eventcellnum,celltype    # 返回所有细胞的突变位点、生长总时间、发生转移时的细胞总数、所有细胞的细胞类型
 
#%% 
b=math.log(2) # 出生率
d=0 # 死亡率
nend=100000# 细胞总数 20微米大小 1厘米可以检测到 那么边长有500个 立体就有125000000个 1000万个 这里用10万个先
s=[0] # 适合度，中性进化
mu=20 # 突变率
ploidy=2 #二倍体
tevent=[8] # 产生亚克隆的时刻，表示早转移和晚转移 #0.1,0.2,0.5,1,2,5,10
read_depth=200 # 测序阅读深度
detectlimit=5/read_depth # 能检测到的位点最低频率
cell,t,ecn,ctp=grow(b=b,d=d,nend=nend,mu=mu,s=s,tevent=tevent)
nummeta=sum(ctp)
#%%
#压扁函数
def flat(l):
    for k in l:
        if not isinstance(k, (list, tuple)):
            yield k
        else:
            yield from flat(k)
            
#计算频数
from collections import Counter

             
def af(cell):    # 获得突变位点的突变频率，序号
    count=Counter(sorted(list(flat(cell))))
    AF=list(count.values())
    return AF,count

def cal_vaf(AF,cellnum,detectlimit=0.01,read_depth=100):    # 模拟测序过程产生vaf
    AF=[x/ploidy for x in AF]
    AF=list(filter(lambda x: x>detectlimit*cellnum,AF))
    depth = list(np.random.poisson(read_depth,len(AF)))
    samp_alleles = list(map(lambda x,y: np.random.binomial(x,y/cellnum),depth,AF))
    vaf=list(map(lambda x,y: x/y,samp_alleles,depth))
    return vaf   
#%%
def aff1(count,AF):  # 获得突变位点的突变频率，有序号 方法一 慢
    AFF=[]
    j=0
    for i in range(0,len(AF)-1):
        if i+1==list(count.keys())[j]:
            AFF.append(list(count.values())[j])
            j=j+1
        else:
            AFF.append(0)
    return AFF
#%%
def aff2(count,AF):  # 获得突变位点的突变频率，有序号 方法二 慢
    AFF=[]
    j=1
    for k,v in count.items():
        while j<k:
            AFF.append(0)
            j=j+1
            if j==k:
                AFF.append(v)
    while j<len(AF):
        AFF.append(0)
    return AFF

#%%
def aff(count,AF_count):
    aa=dict(zip(list(AF_count.keys()),list(0 for i in range(len(AF_count)))))
    AFF={}
    AFF.update(aa)
    AFF.update(count)
    AFF=list(AFF.values())
    return AFF


        
#%% 总的

AF,AF_count=af(cell)
vaf=np.array(AF)/(nend*ploidy)
with open('vaf_all.txt', 'w') as fl:
    for item in vaf:
        fl.write("%s\n" % item)
        
        #%%

## 后面是对细胞进行分类，需要的时候再用上
cell_meta_idx=list(np.where([ctpi==1 for ctpi in ctp])[0])
cell_meta=list(cell[i] for i in cell_meta_idx)
cell_prim_idx=list(np.where([ctpi==0 for ctpi in ctp])[0])
cell_prim=list(cell[i] for i in cell_prim_idx)


#%%

#不考虑深度 用真实值true

AF_meta,AF_meta_count=af(cell_meta)
vaf_meta=np.array(AF_meta)/(nummeta*ploidy)
AF_prim,AF_prim_count=af(cell_prim)
vaf_prim=np.array(AF_prim)/((nend-nummeta)*ploidy)
with open('vaf_meta.txt', 'w') as fl:
    for item in vaf_meta:
        fl.write("%s\n" % item)
with open('vaf_prim.txt', 'w') as fl:
    for item in vaf_prim:
        fl.write("%s\n" % item) 
        
    #%%
# 统计包括零的元素
AFF_prim=aff(AF_prim_count,AF_count)
vaff_prim=np.array(AFF_prim)/((nend-nummeta)*ploidy)
AFF_meta=aff(AF_meta_count,AF_count)
vaff_meta=np.array(AFF_meta)/(nummeta*ploidy)
with open('vaff_meta.txt', 'w') as fl:
    for item in vaff_meta:
        fl.write("%s\n" % item)
with open('vaff_prim.txt', 'w') as fl:
    for item in vaff_prim:
        fl.write("%s\n" % item)
#%%
## 考虑测序深度
vaf=cal_vaf(AF=AF,cellnum=nend,detectlimit=detectlimit,read_depth=read_depth)
with open('vaf_all.txt', 'w') as fl:
    for item in vaf:
        fl.write("%s\n" % item)
vaf_meta=cal_vaf(AF=AF_meta,cellnum=nend,detectlimit=detectlimit,read_depth=read_depth)
with open('vaf_meta.txt', 'w') as fl:
    for item in vaf_meta:
        fl.write("%s\n" % item)

vaf_prim=cal_vaf(AF=AF_prim,cellnum=nend,detectlimit=detectlimit,read_depth=read_depth)
with open('vaf_prim.txt', 'w') as fl:
    for item in vaf_prim:
        fl.write("%s\n" % item) 


