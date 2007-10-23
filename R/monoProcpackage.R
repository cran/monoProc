setOldClass("locfit")
setOldClass("loess")
setClassUnion("fit", c("list","locfit", "loess", "locpoly", "ksmooth"))
setClass("monofit", representation(x="numeric",y="numeric",z="matrix"))
setClass("monoproc",  representation(fit="monofit", fitold="fit", gridsize="numeric", bandwidth="numeric", kernel="character", mono="character", name="character", call="call"))
setClass("monoproc.1d", representation("monoproc"), 
	validity = function(object){
	if(!(length(object@fit@z)==0)) return("monotone fit is not of dimension 1")
	if(object@bandwidth<0) return("bandwidth is not positive")
	if(object@mono!="increasing"&& object@mono!="decreasing") return("Monotonicity constraint is not identifiable")
}
)
setClass("monoproc.2d", representation("monoproc", dir="character"),
validity = function(object){
	#if((length(object@fit)==0)) return("monotone fit is not of dimension 2")
	if(object@bandwidth<0) return("bandwidth is not positive")
	if(length(object@mono)!=2) return("Monotonicity constraint is not identifiable")
	if(object@mono[1]!="increasing"&& object@mono[1]!="decreasing") return("Monotonicity constraint for x is not identifiable")
	if(object@mono[2]!="increasing"&& object@mono[2]!="decreasing") return("Monotonicity constraint for y is not identifiable")
	if(object@dir!="x"&& object@dir!="y"&&object@dir!="xy"&& object@dir!="yx") return("dir is not clear")
})
setClass("monoproclocfit.1d", representation("monoproc.1d"), 
	validity=function(object){
	if(class(object@fitold)!="locfit")cat("this is not an object of class monoproclocfit.1d\n")
	})
setClass("monoproclocfit.2d", representation("monoproc.2d"), 
	validity=function(object){
	if(class(object@fitold)!="locfit")cat("this is not an object of class monoproclocfit.2d\n")
	})

mono.2d<-function(fit, bandwidth,xx, kernel="epanech",dir="xy",mono1, mono2)
{	
	N<-length(fit[[1]])
	n<-length(fit[[2]])
	h<-bandwidth
	dir<-match.arg(dir, c("x", "y", "xy", "yx"))
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	mono2<-match.arg(mono2, c("increasing", "decreasing"))
	if(missing(xx)||class(xx)!="list"){
	 	xx<-list(x=fit[[1]],y=fit[[2]])
		N1=N
		N2=n
	}else{
	N1<-length(xx[[1]])
	N2<-length(xx[[2]])
	}
	if(class(N1)!="integer"||class(N2)!="integer"){
	xx<-list(x=fit[[1]],y=fit[[2]])
		N1=N
		N2=n
	}
	if(dir=="xy"||dir=="yx"){
	if(class(fit[[3]])!="matrix") z=matrix(fit[[3]],ncol=n,nrow=N)  
	else if(class(fit[[3]])=="matrix" && ncol(z)!=n &&nrow(z)!=N){
	 cat("the response does not match to the variables\n")
	z=fit[[3]]
	}}
	else if(dir=="x"){
	if(class(fit[[3]])!="matrix") z=matrix(fit[[3]],ncol=N2,nrow=N)  
	else if(class(fit[[3]])=="matrix" && ncol(z)!=N2 &&nrow(z)!=N){
	 cat("the response does not match to the variables\n") 
	z=fit[[3]]
	}}
	else if(dir=="y"){
	if(class(fit[[3]])!="matrix") z=matrix(fit[[3]],ncol=n,nrow=N1)  
	else if(class(fit[[3]])=="matrix" && ncol(z)!=n &&nrow(z)!=N1){
	 cat("the response does not match to the variables\n") 
	z=fit[[3]]
	}}
	if(dir=="x"){
		zm2<-matrix(0,ncol=N2,nrow=N1)
		for(i in 1:N2)
		{
			fiti<-list(x=fit[[2]],y=z[,i])
			monofit<-mono.1d(fiti,bandwidth=h,xx[[1]],kernel="epanech",mono1=mono1)
			zm2[,i]<-monofit@y
		}	
		fit<-new("monofit",x =xx[[1]], y = xx[[2]],z=zm2)
	}
	else if(dir=="y"){
		zm2<-matrix(0,ncol=N2,nrow=N1)
		for(i in 1:N1)
		{
		fiti<-list(x=fit[[2]],y=z[i,])
		monofit<-mono.1d(fiti,bandwidth=h,xx[[2]],kernel="epanech",mono1=mono2)
		zm2[i,]<-monofit@y
		}	
		fit<-new("monofit",x = xx[[1]], y = xx[[2]],z=zm2)
	}
	else if(dir=="xy"){
		zm<-matrix(0,ncol=n,nrow=N1)
		for(i in 1:n)
		{
			fiti<-list(x=fit[[1]],y=z[,i])
			monofit<-mono.1d(fiti,bandwidth=h,xx[[1]],kernel="epanech",mono1=mono1)
			zm[,i]<-monofit@y
		}	
		zm2<-matrix(0,ncol=N2,nrow=N1)
		for(i in 1:N1)
		{
			fiti2<-list(x=fit[[2]],y=zm[i,])
			monofit<-mono.1d(fiti2,bandwidth=h, xx[[2]],kernel="epanech",mono1=mono2)
			zm2[i,]<-monofit@y
		}	
		fit<-new("monofit",x = xx[[1]], y = xx[[2]],z=zm2)
	}
	else if(dir=="yx"){
		zm<-matrix(0,ncol=N2,nrow=N)
		fiti<-list(x=fit[[2]], y= NULL)
		for(i in 1:N)
		{
			fiti$y = z[i,]
			monofit<-mono.1d(fiti,bandwidth=h, xx[[2]],kernel="epanech",mono1=mono2)
			zm[i,]<-monofit@y
		}	            
		zm2<-matrix(0,ncol=N2,nrow=N1)
		for(i in 1:N2)
		{
			fiti2<-list(x=fit[[1]],y=zm[,i])
			monofit<-mono.1d(fiti2,bandwidth=h, xx[[1]],kernel="epanech",mono1=mono1)
			zm2[,i]<-monofit@y
		}	
		fit<-new("monofit",x = xx[[1]], y = xx[[2]],z=zm2)
	}
	return(fit)

}



"mono.1d"<-function(fit, bandwidth, xx, kernel="epanech", mono1)
{	
	n<-length(fit[[1]])
	h<-bandwidth
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	if(missing(xx)||class(xx)!="numeric") xx<-fit[[1]]
	xl<-(xx-min(fit[[1]]))/(max(fit[[1]])-min(fit[[1]]))
	yl<-fit[[2]]
	N<-length(xl)
	r<-numeric(N)
	monotoneinverse<-function(tn,kernel="epanech")
	{
		xx<-(yl-tn)/h
		g<-rep(1,n)
		if (kernel == "epanech"){
		if(mono1=="increasing") (sum(1/2-3/4*(xx[abs(xx)<1]-1/3*xx[abs(xx)<1]^3))+sum(g[xx<(-1)]))/n
		else if(mono1=="decreasing") (sum(1/2+3/4*(xx[abs(xx)<1]-1/3*xx[abs(xx)<1]^3))+sum(g[xx>1]))/n
		}
	}
	for(l in 1:N)
	{
		tmax  = ifelse(mono1 == "increasing", max(yl), min(yl))
		tmin  = ifelse(mono1 == "increasing", min(yl), max(yl))
		wertmin<-monotoneinverse(tmin)
		wertmax<-monotoneinverse(tmax)
		if(wertmax==wertmin) r[l]<-tmin
		wertt<-xl[l]
		for(i in 1:150)
		{
			tneu<-tmin+(wertt-wertmin)*(tmax-tmin)/(wertmax-wertmin)
			wertneu<-monotoneinverse(tneu)
			if(abs(wertt-wertneu)<0.0001) break;
			if(tneu<(min(yl)-0.1)) break;
			if(tneu>(max(yl)+0.1)) break;
			if(wertneu >wertt)
			{
				tmax<- tneu
				wertmax<-wertneu
				tmin <- tmin
				wertmin <-wertmin
			} else 
			{
				tmin <- tneu
				wertmin <- wertneu
				wertmax <-wertmax
				tmax<-tmax
			}
		}
	
		r[l]<-(tmin + (wertt-wertmin)*(tmax - tmin)/(wertmax- wertmin))
	}
	fit<-new("monofit",x = xx, y = r)
	return(fit)
}



setGeneric("monoproc", function(fit, bandwidth, xx, dir, ...)standardGeneric("monoproc"))
setGeneric("cv", function(fit,...) standardGeneric("cv"))

setMethod("monoproc", signature(fit="locfit", bandwidth="numeric", xx="missing", dir="missing"), function(fit, bandwidth, kernel="epanech", mono1="increasing", gridsize=40){
	call<-match.call(call=sys.call(sys.parent()))
	h<-bandwidth
	ufit<-lfmarg(fit)
	if(length(ufit)!=1)stop("dir is not specified\n")
	n<-gridsize
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	newdata<-(max(ufit[[1]])-min(ufit[[1]]))*seq(1:n)/(n+1)+min(ufit[[1]])
	z=predict(fit,newdata)
	fiti<-list(x=newdata,y=z)
	mopro<-mono.1d(fiti,bandwidth=h,xx=newdata,kernel=kernel,mono1=mono1) 
	monoproc<-new("monoproclocfit.1d", fit=mopro,fitold=fit, gridsize=n, bandwidth=h,kernel=kernel, mono=mono1, name=c(fit$vnames, fit$yname), call=call)
return(monoproc)
}
)
setMethod("monoproc", signature(fit="locfit", bandwidth="numeric", xx="numeric", dir="missing"), function(fit, bandwidth,xx, kernel="epanech", mono1="increasing", gridsize=40){
	call<-match.call(call=sys.call(sys.parent()))
	h<-bandwidth
	ufit<-lfmarg(fit)
	if(length(ufit)!=1)stop("dir is not specified!\n")
	n<-gridsize
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	newdata<-(max(ufit[[1]])-min(ufit[[1]]))*seq(1:n)/(n+1)+min(ufit[[1]])
	z=predict(fit,newdata)
			fiti<-list(x=newdata,y=z)
			mopro<-mono.1d(fiti,bandwidth=h,xx=xx,kernel=kernel,mono1=mono1) 
			monoproc<-new("monoproclocfit.1d", fit=mopro,fitold=fit, gridsize=n, bandwidth=h,kernel=kernel, mono=mono1, name=c(fit$vnames, fit$yname), call=call)
return(monoproc)
}
)

setMethod("monoproc", signature(fit="fit", bandwidth="missing", xx="ANY", dir="ANY"), function(fit, bandwidth,xx, dir,...){
cat("the bandwidth for the monotonizing has to be specified!\n")
}
)

setMethod("monoproc", signature(fit="locfit", bandwidth="numeric", xx="missing", dir="character"), function(fit, bandwidth,xx, dir, kernel="epanech", mono1="increasing", mono2="increasing", gridsize=40){
	call<-match.call(call=sys.call(sys.parent()))
	h<-bandwidth
	ufit<-lfmarg(fit)
	if(length(ufit)!=2) stop("the locfit object has not the right dimension\n
				if the independent variable x of this locfit object is only of length 1, dir may not be used")
	n<-gridsize
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	mono2<-match.arg(mono2, c("increasing", "decreasing"))
	dir<-match.arg(dir, c("x", "y", "xy", "yx"))
	newdata<-list(x=(max(ufit[[1]])-min(ufit[[1]]))*seq(1:n)/(n+1)+min(ufit[[1]]), y=(max(ufit[[2]])-min(ufit[[2]]))*seq(1:n)/(n+1)+min(ufit[[2]]))
	z=predict(fit,newdata)
	fiti<-list(x=newdata$x,y=newdata$y,z=z)
	mopro<-mono.2d(fiti,bandwidth=h,xx=list(x=newdata$x, y=newdata$y), kernel=kernel, dir=dir, mono1=mono1, mono2=mono2)
	monoproc<-new("monoproclocfit.2d", fit=mopro,fitold=fit, gridsize=n, bandwidth=h, kernel=kernel, mono=c(mono1,mono2), dir=dir, name=c(fit$vnames[1], fit$vnames[2], fit$yname), call=call)
return(monoproc)
}
)
setMethod("monoproc", signature(fit="locfit", bandwidth="numeric", xx="list", dir="character"), function(fit, bandwidth,xx, dir, kernel="epanech", mono1="increasing", mono2="increasing", gridsize=40){
	call<-match.call(call=sys.call(sys.parent()))
	h<-bandwidth
	ufit<-lfmarg(fit)
	if(length(ufit)!=2) stop("the locfit object has not the right dimension\n
				if the independent variable x of this locfit object is only of length 1, dir may not be used")
	n<-gridsize
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	mono2<-match.arg(mono2, c("increasing", "decreasing"))
	dir<-match.arg(dir, c("x", "y", "xy", "yx"))
	if(dir=="x"){
		newdata<-list(x=(max(ufit[[1]])-min(ufit[[1]]))*seq(1:n)/(n+1)+min(ufit[[1]]), y=xx[[2]])
	} else if(dir=="y"){
		newdata<-list(x=xx[[1]], y=(max(ufit[[2]])-min(ufit[[2]]))*seq(1:n)/(n+1)+min(ufit[[2]]))
	} else { 
		newdata<-list(x=(max(ufit[[1]])-min(ufit[[1]]))*seq(1:n)/(n+1)+min(ufit[[1]]),y=(max(ufit[[2]])-min(ufit[[2]]))*seq(1:n)/(n+1)+min(ufit[[2]]))
	}
	z=predict(fit,newdata)
	fiti<-list(x=newdata$x,y=newdata$y,z=z)
	mopro<-mono.2d(fiti,bandwidth=h,xx=xx,kernel=kernel,dir=dir,mono1=mono1, mono2=mono2)
	monoproc<-new("monoproclocfit.2d", fit=mopro,fitold=fit, gridsize=n, bandwidth=h,kernel=kernel, mono=c(mono1,mono2), dir=dir, name=c(fit$vnames[1], fit$vnames[2], fit$yname), call=call)
return(monoproc)
}
)



setMethod("monoproc", signature(fit="loess", bandwidth="numeric", xx="missing", dir="missing"), function(fit, bandwidth, kernel="epanech", mono1="increasing", gridsize)
{
	h<-bandwidth
	n<-gridsize
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	if(ncol(fit$x)!=1) stop("monoproc currently only implemented for loess objects if x is of length 1\n")
	newdata<-(max(fit$x)-min(fit$x))*seq(1:n)/(n+1)+min(fit$x)
	z=predict(fit, newdata)
	fiti<-list(x=newdata,y=z)
	if(ncol(fit$x)==1) mopro<-mono.1d(fiti,h,xx=newdata,kernel=kernel,mono1=mono1)
	monoproc<-new("monoproc.1d", fit=mopro,fitold=fit,gridsize=n, bandwidth=h,kernel=kernel, mono=mono1, name=c(as.character(fit$call[[2]][3]), as.character(fit$call[[2]][2])),call=match.call(call=sys.call(sys.parent())))
	return(monoproc)
}
)
setMethod("monoproc", signature(fit="loess", bandwidth="numeric", xx="numeric", dir="missing"), function(fit, bandwidth,xx, kernel="epanech", mono1="increasing", gridsize)
{
	h<-bandwidth
	n<-gridsize
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	if(ncol(fit$x)!=1)stop("monoproc currently only implemented for loess objects if x is of length 1\n")
	newdata<-(max(fit$x)-min(fit$x))*seq(1:n)/(n+1)+min(fit$x)
	z=predict(fit, newdata)
	fiti<-list(x=newdata,y=z)
	if(ncol(fit$x)==1) mopro<-mono.1d(fiti,h,xx=xx,kernel=kernel,mono1=mono1)
	monoproc<-new("monoproc.1d", fit=mopro,fitold=fit,gridsize=n, bandwidth=h,kernel=kernel, mono=mono1, name=c(as.character(fit$call[[2]][3]), as.character(fit$call[[2]][2])),call=match.call(call=sys.call(sys.parent())))
	return(monoproc)
}
)


setMethod("monoproc", signature(fit="locpoly", bandwidth="numeric", xx="missing", dir="missing"), function(fit, bandwidth, kernel="epanech",mono1="increasing", gridsize)
{
	h<-bandwidth
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	if(missing(gridsize)) {
	n<-length(fit@x)
	fitnew<-list(x=fit@x, y=fit@y)
	}else {
	n<-gridsize
	fit@call$gridsize<-n
	pfitnew<-eval(fit@call)
	fitnew<-list(x=pfitnew@x, y=pfitnew@y)}
	mopro<-mono.1d(fitnew,h,xx=fit@x,kernel="epanech",mono1=mono1)
	monoproc<-new("monoproc.1d", fit=mopro,fitold=fit,gridsize=n, bandwidth=h,kernel=kernel, mono=mono1,name=c(as.character(fit@call[[2]]), as.character(fit@call[[3]])),call=match.call(call=sys.call(sys.parent())))
	return(monoproc)
})
setMethod("monoproc", signature(fit="locpoly", bandwidth="numeric", xx="numeric", dir="missing"), function(fit, bandwidth,xx, kernel="epanech",mono1="increasing", gridsize)
{
	h<-bandwidth
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	if(missing(gridsize)) {
	n<-length(fit@x)
	fitnew<-list(x=fit@x, y=fit@y)
	}else {
	n<-gridsize
	fit@call$gridsize<-n
	pfitnew<-eval(fit@call)
	fitnew<-list(x=pfitnew@x, y=pfitnew@y)}
	mopro<-mono.1d(fitnew,h,xx,kernel="epanech",mono1=mono1)
	monoproc<-new("monoproc.1d", fit=mopro,fitold=fit,gridsize=n, bandwidth=h,kernel=kernel, mono=mono1,name=c(as.character(fit@call[[2]]), as.character(fit@call[[3]])),call=match.call(call=sys.call(sys.parent())))
	return(monoproc)
})

setMethod("monoproc", signature(fit="ksmooth", bandwidth="numeric", xx="missing", dir="missing"), function(fit, bandwidth, kernel="epanech", mono1="increasing",gridsize)
{
	h<-bandwidth
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	if(missing(gridsize)) {
	n<-length(fit@x)
	fitnew<-list(x=fit@x, y=fit@y)
	}else {
	n<-gridsize
	fit@call$n.points<-gridsize
	pfitnew<-eval(fit@call)
	fitnew<-list(x=pfitnew@x, y=pfitnew@y)
	}
	mopro<-mono.1d(fitnew,h,xx=fit@x,kernel="epanech",mono1=mono1)
	monoproc<-new("monoproc.1d", fit=mopro,fitold=fit,gridsize=n, bandwidth=h,kernel=kernel, 			mono=mono1,name=c(as.character(fit@call[[2]]), as.character(fit@call[[3]])),call=match.call(call=sys.call(sys.parent())))
	return(monoproc)
})
setMethod("monoproc", signature(fit="ksmooth", bandwidth="numeric", xx="numeric", dir="missing"), function(fit, bandwidth,xx, kernel="epanech", mono1="increasing",gridsize)
{

	h<-bandwidth
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	if(missing(gridsize)) {
	n<-length(fit@x)
	fitnew<-list(x=fit@x, y=fit@y)
	}else {
	n<-gridsize
	fit@call$n.points<-gridsize
	pfitnew<-eval(fit@call)
	fitnew<-list(x=pfitnew@x, y=pfitnew@y)
	}
	mopro<-mono.1d(fitnew,h,xx,kernel="epanech",mono1=mono1)
	monoproc<-new("monoproc.1d", fit=mopro,fitold=fit,gridsize=n, bandwidth=h,kernel=kernel, 			mono=mono1,name=c(as.character(fit@call[[2]]), as.character(fit@call[[3]])),call=match.call(call=sys.call(sys.parent())))
	return(monoproc)
})





setMethod("monoproc", signature(fit="list", bandwidth="numeric", xx="missing", dir="missing"), function(fit, bandwidth, kernel="epanech", mono1="increasing")
{
	if(length(fit)!=2)stop("This is not a fit of length 2! \n If it is a fit of length 3, dir has to be specified!\n")
	h<-bandwidth
	n=length(fit[[1]])
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	name=c(names(fit[1]), names(fit[2]))
	if(is.null(name)) name=c("x","y")
	mopro<-mono.1d(fit,h,xx=fit[[1]],kernel="epanech",mono1=mono1)
	monoproc<-new("monoproc.1d", fit=mopro,fitold=fit,gridsize=n, bandwidth=h,kernel=kernel,mono=mono1,name=name, call=match.call(call=sys.call(sys.parent())))
	return(monoproc)
})
setMethod("monoproc", signature(fit="list", bandwidth="numeric", xx="numeric", dir="missing"), function(fit, bandwidth,xx, kernel="epanech", mono1="increasing")
{
	if(length(fit)!=2)stop("This is not a fit of length 2!\n If it is a fit of length 3, dir has to be specified!\n")
	h<-bandwidth
	n=length(fit[[1]])
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	name=c(names(fit[1]),names(fit[2]))
	if(is.null(name)) name=c("x","y")
	mopro<-mono.1d(fit,h,xx,kernel="epanech",mono1=mono1)
	monoproc<-new("monoproc.1d", fit=mopro,fitold=fit,gridsize=n, bandwidth=h,kernel=kernel, mono=mono1,name=name,call=match.call())
	return(monoproc)
})
setMethod("monoproc", signature(fit="list", bandwidth="numeric", xx="missing", dir="character"), function(fit, bandwidth,dir, kernel="epanech", mono1="increasing",mono2="increasing")
{
	if(length(fit)!=3) stop("This is not a fit of length 3!\n")
	h<-bandwidth
	n=length(fit[[1]])
	N=length(fit[[2]])
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	mono2<-match.arg(mono2, c("increasing", "decreasing"))
	name= c(names(fit[1]), names(fit[2]), names(fit[3]))
	if(is.null(name)) name<-c("x","y","z")
	dir<-match.arg(dir, c("x", "y", "xy", "yx"))
	mopro<-mono.2d(fit,h,xx=list(x=fit[[1]],y=fit[[2]]),kernel=kernel,dir=dir,mono1=mono1, mono2=mono2)
	monoproc<-new("monoproc.2d", fit=mopro,fitold=fit,gridsize=c(n,N), bandwidth=h,kernel=kernel, mono=c(mono1,mono2), dir=dir, name=name, call=match.call())
	return(monoproc)
})
setMethod("monoproc", signature(fit="list", bandwidth="numeric", xx="list", dir="character"), function(fit, bandwidth,xx, dir, kernel="epanech", mono1="increasing",mono2="increasing")
{
	if(length(fit)!=3)stop("This is not a fit of length 3!\n")
	h<-bandwidth
	n=length(fit[[1]])
	N=length(fit[[2]])
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	mono2<-match.arg(mono2, c("increasing", "decreasing"))
	name= c(names(fit[1]), names(fit[2]), names(fit[3]))
	if(is.null(name)) name<-c("x","y","z")
	dir<-match.arg(dir, c("x", "y", "xy", "yx"))
	mopro<-mono.2d(fit,h,xx=xx,kernel=kernel,dir=dir,mono1=mono1, mono2=mono2)
	monoproc<-new("monoproc.2d", fit=mopro,fitold=fit,gridsize=c(n,N), bandwidth=h,kernel=kernel, mono=c(mono1,mono2), dir=dir, name=name, call=match.call(call=sys.call(sys.parent())))
	return(monoproc)
})

setMethod("plot", signature(x="monofit",y="missing"), function(x,y,...){
plot(x@x,x@y,...)
})

setMethod("lines", signature(x="monofit"), function(x,...){
lines(x@x,x@y,...)
})

setMethod("plot", signature(x="monoproc.1d",y="missing"), function(x,xlab=x@name[1],ylab=x@name[2],...){
plot(x@fit@x,x@fit@y,type="l",xlab=x@name[1],ylab=x@name[2],...)
})

setMethod("lines", signature(x="monoproc.1d"), function(x,...){
lines(x@fit@x,x@fit@y,...)
})

setMethod("plot",signature(x="monoproc.2d", y="missing"), function(x,type="persp",xlab=x@name[1],ylab=x@name[2],zlab=x@name[3],...){ 
if(type=="persp") persp(x@fit@x,x@fit@y,x@fit@z,xlab=x@name[1],ylab=x@name[2],zlab=x@name[3],...)
if(type=="contour") contour(x@fit@x,x@fit@y,x@fit@z,xlab=x@name[1],ylab=x@name[2],...)
if(type=="image") image(x@fit@x,x@fit@y,x@fit@z,xlab=x@name[1],ylab=x@name[2],...)
})


setMethod("summary",signature(object="monoproc.1d"), function(object,...)
{
if(!(class(object@fitold)=="list")){
	cat("Old Model:\n")
	print(object@fitold)
}
cat("Smooth monotone", object@mono, "fit\n")
cat("Smoothing parameter:	", object@bandwidth,"\n")
cat("Kernel:	", object@kernel, "\n") 
cat("N:	",object@gridsize, "\n")
cat("Evaluated points:	", length(object@fit@x)*length(object@fit@y),"\n")
cat("fitted values:\n")
object@fit
}
)

setMethod("summary",signature(object="monoproc.2d"), function(object,...)
{
if(!(class(object@fitold)=="list")){
	cat("Old Model:\n")
	print(object@fitold)
}
if(object@dir=="xy"||object@dir=="yx") {
cat("Smooth monotone", object@mono[1], "with respect to x and monotone", object@mono[2], "with respect to y fit\n")
cat("Order of monotonization: ", object@dir, "\n")}
else if(object@dir=="x") cat("Smooth monotone", object@mono[1], "with respect to x fit\n")
else if(object@dir=="y") cat("Smooth monotone", object@mono[2], "with respect to y fit\n")
cat("Smoothing parameter:	", object@bandwidth,"\n")
cat("Kernel:	", object@kernel, "\n") 
cat("N:	",object@gridsize, "\n")
cat("Evaluated points:	", length(object@fit@x)*length(object@fit@y),"\n")
cat("fitted values:\n")
object@fit
}
)


setMethod("print",signature(x="monoproc.1d"), function(x,...)
{
if(!(class(x@fitold)=="list")){
	cat("Old Model:\n")
	print(x@fitold)

}

cat("Smooth monotone", x@mono, "fit\n")
cat("Smoothing parameter:	", x@bandwidth,"\n")
cat("Kernel:	", x@kernel, "\n") 
cat("N:	",x@gridsize, "\n")
cat("Evaluated points:	", length(x@fit@x),"\n")
cat("the first fitted value:\n")
print(list(x=x@fit@x[1],y=x@fit@y[1]))
}
)

setMethod("show",signature(object="monoproc.1d"), function(object)
{
if(!(class(object@fitold)=="list")){
	cat("Old Model:\n")
	print(object@fitold)

}

cat("Smooth monotone", object@mono, "fit\n")
cat("Smoothing parameter:	", object@bandwidth,"\n")
cat("Kernel:	", object@kernel, "\n") 
cat("N:	",object@gridsize, "\n")
cat("Evaluated points:	", length(object@fit@x),"\n")
cat("the first fitted value:\n")
print(list(x=object@fit@x[1],y=object@fit@y[1]))
}
)

setMethod("print",signature(x="monoproc.2d"), function(x,...)
{
if(!(class(x@fitold)=="list")){
	cat("Old Model:\n")
	print(x@fitold)

}

if(x@dir=="xy"||x@dir=="yx") {
cat("Smooth monotone", x@mono[1], "with respect to x and monotone", x@mono[2], "with respect to y fit\n")
cat("Order of monotonization: ", x@dir, "\n")}
else if(x@dir=="x") cat("Smooth monotone", x@mono[1], "with respect to x fit\n")
else if(x@dir=="y") cat("Smooth monotone", x@mono[2], "with respect to y fit\n")
cat("Smoothing parameter:	", x@bandwidth,"\n")
cat("Kernel:	", x@kernel, "\n")
cat("N:	",x@gridsize, "\n")
cat("Evaluated points:	", length(x@fit@x)*length(x@fit@y),"\n")
cat("the first fitted value:\n")
print(list(x=x@fit@x[1],y=x@fit@y[1],z=x@fit@z[1,1]))
}
)

setMethod("show",signature(object="monoproc.2d"), function(object)
{
if(!(class(object@fitold)=="list")){
	cat("Old Model:\n")
	print(object@fitold)

}

if(object@dir=="xy"||object@dir=="yx") {
cat("Smooth monotone", object@mono[1], "with respect to x and monotone", object@mono[2], "with respect to y fit\n")
cat("Order of monotonization: ", object@dir, "\n")}
else if(object@dir=="x") cat("Smooth monotone", object@mono[1], "with respect to x fit\n")
else if(object@dir=="y") cat("Smooth monotone", object@mono[2], "with respect to y fit\n")
cat("Smoothing parameter:	", object@bandwidth,"\n")
cat("Kernel:	", object@kernel, "\n")
cat("N:	",object@gridsize, "\n")
cat("Evaluated points:	", length(object@fit@x)*length(object@fit@y),"\n")
cat("the first fitted value:\n")
print(list(x=object@fit@x[1],y=object@fit@y[1],z=object@fit@z[1,1]))
}
)

setMethod("cv", signature(fit="monoproc.1d"), function(fit){	
	fitold<-fit@fitold
       data <-{ if (is.null(fitold$call$data)){ 
           sys.frame(sys.parent())}        else eval(fitold$call$data)}
 	m<-locfit.matrix(fitold, data = data)
	n<-length(m$y)
	cv<-matrix(0,nrow=n, ncol=2)
	cv[,2]<-fitted.locfit(fitold,cv=TRUE)
	for(i in 1:n){
	fitold$call[[2]]<-as.matrix(m$x[-i,])
	fitold$call[[3]]<-as.numeric(m$y[-i])
	fitcv<-eval(fitold$call)
	t<-monoproc(fitcv,bandwidth=fit@bandwidth,  xx=m$x[i,1], kernel=fit@kernel,mono1=fit@mono,gridsize=fit@gridsize)@fit
	cv[i,1]<-t@y
	}
	return(cv)
	}
)

setMethod("cv", signature(fit="monoproclocfit.2d"), function(fit,...){
	
	fitold<-fit@fitold
        	data <-{ if (is.null(fitold$call$data)){ 
          	sys.frame(sys.parent())}
       	else eval(fitold$call$data)}
  	m<-locfit.matrix(fitold, data = data)
	n<-length(m$y)
	cv<-matrix(0,nrow=n, ncol=2)
	cv[,2] <- fitted.locfit(fitold,cv=TRUE)
	for(i in 1:n){
		fitold$call[[2]] <- as.matrix(m$x[-i,])
		fitold$call[[3]] <- as.numeric(m$y[-i])
		fitcv <- eval(fitold$call)
		t <- monoproc(fitcv, bandwidth=fit@bandwidth,
                               xx=list(x=m$x[i,1],y=m$x[i,2]), dir=fit@dir,
                               kernel=fit@kernel,mono1=fit@mono[[1]], mono2=fit@mono[[2]], gridsize=fit@gridsize)@fit
		cv[i,1] <- t@z
	}
	return(cv)
	}
)




