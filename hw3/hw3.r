data=c(404,96,301,199,348,152,235,265,287,213,183,317,251,249,102,398,148,352,20,480)

nll.condition=function(par,y) { 
    p=1:4
    d=par[1]
    g=par[2]
    p[1]=d+(1-d)*g
    p[2]=1-p[1]
    p[3]=g
    p[4]=1-g
    return(-sum(y*log(p)))
}

nll.binomial=function(par10,dat) {
    model=1:10
    result=0.0
    par=1:2
    params=NULL
    for (i in 0:4) {
        start2=i*2+1
        end2=start2+1
        start4=i*4+1
        end4=start4+3
        model=optim(par=par10[start2:end2], fn=nll.condition, y=dat[start4:end4])
        result=-lchoose(dat[start4]+dat[start4+1], dat[start4])-lchoose(dat[start4+2]+dat[start4+3], dat[start4+3])+result+model$value
        params=c(params,model$par)
    }
    return(result)
}

nll.highThreshold.given.d=function(g,y,d) {
    p=1:4
    p[1]=d+(1-d)*g 
    p[2]=1-p[1]
    p[3]=g 
    p[4]=1-g
    return(-sum(y*log(p)))
}
nll.highThreshold=function(d, dat) {
    min=0.0
    for (i in 0:4) {
        start2=i*2+1
        end2=start2+1
        start4=i*4+1
        end4=start4+3
        min=min+optimize(nll.highThreshold.given.d, interval=c(0,1), y=dat[start4:end4], d=d)$objective
        min=min-lchoose(dat[start4]+dat[start4+1], dat[start4])-lchoose(dat[start4+2]+dat[start4+3], dat[start4+3])
    }
    return(min)
}
nll.gerneralHighThreshold.1=function(g,y,det) {
    ds=det[1]
    dn=det[2]
    p=1:4
    p[1]=ds+(1-ds)*g 
    p[2]=1-p[1]
    p[3]=(1-dn)*g 
    p[4]=1-p[3]
    return(-sum(y*log(p)))
}

nll.gerneralHighThreshold=function(zdet,dat) {
    det=plogis(zdet)
    result=0.0
    for (i in 0:4) {
        start4=i*4+1
        end4=start4+3
        result=result+optimize(nll.gerneralHighThreshold.1, interval=c(0,1), y=dat[start4:end4], det=det)$objective
        result=result-lchoose(dat[start4]+dat[start4+1], dat[start4])-lchoose(dat[start4+2]+dat[start4+3], dat[start4+3])
    }
    return(result)
}
ds.amount.of.bias=function(b,p) {
    ifelse(b>0, (1-p)*b, p*b)
}

nll.lowThreshold.1=function(b,y,det) {
    ds=det[1]
    dn=det[2]
    p=1:4
    p[1]=ds.amount.of.bias(b,ds)
    p[2]=1-p[1]
    p[3]=dn+ds.amount.of.bias(b,dn)
    p[4]=1-p[3]
    return(-sum(y*log(p)))
}
 
nll.lowThreshold=function(det,dat) {
    result=0.0
    for (i in 0:4) {
        start4=i*4+1
        end4=start4+3
        result=result+optimize(nll.lowThreshold.1, interval=c(0,1), y=dat[start4:end4], det=det)$objective
        result=result-lchoose(dat[start4]+dat[start4+1], dat[start4])-lchoose(dat[start4+2]+dat[start4+3], dat[start4+3])
    }
    return(result)
}

nll.doubleHighThreshold.1=function(g,y,d) {
    p=1:4
    p[1]=d+(1-d)*g 
    p[2]=1-p[1]
    p[3]=(1-d)*g 
    p[4]=1-p[3]
    return(-sum(y*log(p)))
}
nll.doubleHighThreshold=function(d,dat) {
    min=0.0
    for (i in 0:4) {
        start2=i*2+1
        end2=start2+1
        start4=i*4+1
        end4=start4+3
        min=min+optimize(nll.doubleHighThreshold.1, interval=c(0,1), y=dat[start4:end4], d=d)$objective
        min=min-lchoose(dat[start4]+dat[start4+1], dat[start4])-lchoose(dat[start4+2]+dat[start4+3], dat[start4+3])
    }
    return(min)
}

nll.freeVariance.1=function(c,y,par) {
    d=par[1]
    sigma=par[2]
    p=1:4
    c=qlogis(c)
    p[1]=1-pnorm(c,d,sigma)
    p[2]=1-p[1]
    p[3]=1-pnorm(c)
    p[4]=1-p[3]
    return(-sum(y*log(p)))
}

nll.freeVariance=function(par,dat) {
    result=0.0
    for (i in 0:4) {
        start4=i*4+1
        end4=start4+3
        result=result+optimize(nll.freeVariance.1, interval=c(0,1), y=dat[start4:end4], par=par)$objective
        result=result-lchoose(dat[start4]+dat[start4+1], dat[start4])-lchoose(dat[start4+2]+dat[start4+3], dat[start4+3])
    }
    return(result)
}

nll.equalVariance.1=function(c,y,par) {
    d=par
    sigma=1
    p=1:4
    c=qlogis(c)
    p[1]=1-pnorm(c,d,sigma)
    p[2]=1-p[1]
    p[3]=1-pnorm(c)
    p[4]=1-p[3]
    return(-sum(y*log(p)))
}

nll.equalVariance=function(par,dat) {
    result=0.0
    for (i in 0:4) {
        start4=i*4+1
        end4=start4+3
        result=result+optimize(nll.equalVariance.1, interval=c(0,1), y=dat[start4:end4], par=par)$objective
        result=result-lchoose(dat[start4]+dat[start4+1], dat[start4])-lchoose(dat[start4+2]+dat[start4+3], dat[start4+3])
    }
    return(result)
}

suppressWarnings({

par10=rep(0.5,10)
est.binominal=nll.binomial(par10, dat=data)
est.highThreshold=optimize(nll.highThreshold, interval=c(0,1),dat=data)
zdet=rep(0,2)
est.gernalHighThreshold=optim(zdet, nll.gerneralHighThreshold, dat=data)
par2=rep(0.5,2)
est.freeVariance=optim(par2, nll.freeVariance, dat=data)
par2=rep(0.5,2)
est.lowThreshold=optim(par2, nll.lowThreshold, dat=data)
est.doubleHighThreshold=optimize(nll.doubleHighThreshold, interval=c(0,1),dat=data)
est.equalVariance=optimize(nll.equalVariance, interval=c(0,1), dat=data)

})
logLike=1:7
logLike[1]=est.binominal
logLike[2]=est.gernalHighThreshold$value
logLike[3]=est.highThreshold$objective
logLike[4]=est.doubleHighThreshold$objective
logLike[5]=est.lowThreshold$value
logLike[6]=est.freeVariance$value
logLike[7]=est.equalVariance$objective

models=1:7
models[1]="binominal"
models[2]="high threshold"
models[3]="general threshold"
models[4]="free variance"
models[5]="low threshold"
models[6]="double high threshold"
models[7]="equal variance"

G_squar=1:6
cal_G2=function(x, y) {
    return(2*(logLike[x]-logLike[y]))
}
G_squar[1]=cal_G2(2,1)
G_squar[2]=cal_G2(3,2)
G_squar[3]=cal_G2(4,1)
G_squar[4]=cal_G2(5,2)
G_squar[5]=cal_G2(6,1)
G_squar[6]=cal_G2(7,6)

params=1:7
params[1]=10
params[2]=7
params[3]=6
params[4]=6
params[5]=7
params[6]=7
params[8]=6
AIC=1:7
cal_AIC=function(index) {
    return(2*logLike[index]+2*params[index])
}
for (i in 1:7) {
    AIC[i]=cal_AIC(i)
}
print(models)
print("negative log likelihood")
print(logLike)
print("G squar")
print(G_squar)
print("AIC")
print(AIC)
