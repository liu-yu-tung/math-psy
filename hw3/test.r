suppressWarnings({
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

# same as modeling book
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

nll.freeVariance.1=function(c,y,par) {
    d=par[1]
    simga=par[2]
    p=1:4
    c=qlogis(c)
    p[1]=1-pnorm(c,d,simga)
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



par10=rep(0.5,10)
zdet=rep(0,2)
par2=rep(0.5,3)
est.highThreshold=optimize(nll.highThreshold, interval=c(0,1),dat=data)
est.binominal=nll.binomial(par10, dat=data)
est.gernalHighThreshold=optim(zdet, nll.gerneralHighThreshold, dat=data)
est.freeVariance=optim(par2, nll.freeVariance, dat=data)
print(nll.freeVariance(par=rep(0.8,2),dat=data))
print(est.binominal)
print(est.highThreshold$objective)
print(est.gernalHighThreshold$value)
print(est.freeVariance$value)
})