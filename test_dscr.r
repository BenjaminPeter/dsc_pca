library(dscr)
library(MASS)

source("pcasel_datamakers.r")
source("pcasel_methods.r")
source("pcasel_scores.r")


source("dscr_groups.r")

dsc_pcasel = new.dsc("dsc_pcasel", "dsc_pcasel")


if(F){
addScenario(dsc_pcasel, name="sel1", fn=datamaker,
            args=list(fun=datamaker.mvngenotypes,
                      n.samples=100,
                      n.neutral.snps=9500,
                      n.selected.snps=0500,
                      geno=T)
            , seed=1:2)

addScenario(dsc_pcasel, name="selcos",
            fn=datamaker,
            args=list(fun=datamaker.discrete.cosine,
                      n.samples.per.dim=30,
                      n.neutral.snps=9500,
                      n.selected.snps=0500,
                      geno=T),
            seed=1:2)
addScenario(dsc_pcasel, name="selcos2",
            fn=datamaker, 
            args=list(fun=datamaker.discrete.cosine,
                      n.samples.per.dim=30,
                      n.neutral.snps=9500,
                      n.selected.snps=0500,
                      geno=F),
            seed=1:2)
}

n.samples.per.dim=c(10,20,30)
y.scale <- c(1,2,5)
x.scale <- 1
sel.strength <- c(0,1,5,10)
sel.angle <- seq(0,pi/2, length.out=5)
spread <- c(0.1,0.2,1)
flexargs <- expand.grid(n.samples.per.dim=n.samples.per.dim, 
                        x.scale=x.scale, y.scale=y.scale,
                        sel.strength=sel.strength, 
                        sel.angle=sel.angle,
                        spread=spread)
addScenarioGroup(dsc_pcasel, name="selcos", 
                 fn=datamaker,
                 seed=c(100,200,300),
                 args=list(fun=datamaker.discrete.cosine2,
                           n.neutral.snps=1900,
                           n.selected.snps=100,
                           geno=T),
                 flexible.args=flexargs)

addMethod(dsc_pcasel, name="df.rho1", 
          fn=method.dufouret.rho1)
addMethod(dsc_pcasel, name="df.rho2",
          fn=method.dufouret.rho2)
addMethod(dsc_pcasel, name="df.h", 
          fn=method.dufouret.h,
          args=list(K=2))
addMethod(dsc_pcasel, name="df.hprime",
          args=list(K=2),
          fn=method.dufouret.hprime)

addMethod(dsc_pcasel, name="lfa3",
          args=list(K=3),
          fn=method.lfa)
addMethod(dsc_pcasel, name="fst",
          args=list(K=3),
          fn=method.fst)
#addMethod(dsc_pcasel, name="lfa5",
#          args=list(K=5),
#          fn=method.lfa)
addMethod(dsc_pcasel, name="gwish",
          args=NULL,
          fn=method.gwish)
addMethod(dsc_pcasel, name="norm",
          args=NULL,
          fn=method.normal)


addScore(dsc_pcasel, name="top100", fn=score.top100)
#addScore(dsc_pcasel, name="top50", fn=score.top50)
#addScore(dsc_pcasel, name="top500", fn=score.top500)
