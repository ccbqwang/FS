
FD <- read.table("F.txt")
SD <- read.table("S.txt")



source("inter_sup.R")
source("inter_plot.R")

### data process
dat_i <- dat_index(f=FD,s=SD)

### allo fit
res1_i <- allo_fit(x=dat_i$x,y=dat_i$y,nt=dat_i$times)

Figure11(ret=res1_i,ii=c(1,12,23,44),stage=dat_i$stage,oii=dat_i$orderi)
Figure2(ret=res1_i,ii=c(1,12,23,44),stage=dat_i$stage,oii=dat_i$orderi)

Figure21(ret=res1_i,ii=c(1,12,23,44),stage=dat_i$stage,oii=dat_i$orderi)



Figure1_rec(ret=res1_i,ii=c(1,12,23,44),stage=dat_i$stage,oii=dat_i$orderi)


Figure2_rec(ret=res1_i,ii=c(1,12,23,44),stage=dat_i$stage,oii=dat_i$orderi)


###ode fit  
  res2_i <- ode_fit(dat=res1_i)

###ode_effect_type
ode_eff <- effect_i(dat=res1_i,ret=res2_i,ntl=50)

Figure2_c11(effect=ode_eff)

GO <- GO_enrich_type(effect=ode_eff)
a11 <- read.csv("GO_KEGG/C1.csv")
a22 <- read.csv("GO_KEGG/C2(1).csv")
a33 <- read.csv("GO_KEGG/C3(1).csv")
a44 <- read.csv("GO_KEGG/C4(1).csv")
a55 <- read.csv("GO_KEGG/C5(1).csv")

gosig <- list()
gosig[[1]] <- a11;gosig[[2]] <- a22;gosig[[3]] <- a33;gosig[[4]] <- a44;gosig[[5]] <- a55

Figure2_c2(GO1=gosig,GO2=GO)




gene <- read.csv("../Pro1_6_1/Gene_l1.csv")

match(gene[,1],ode_eff$gname)

type <- ode_eff$type[match(gene[,1],ode_eff$gname)]

ngene1 <- cbind(gene,type=type)
write.csv(ngene1,file="nege1_l.csv")
