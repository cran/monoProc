setClass("ksmooth", representation(x="numeric",y="numeric", call="call"))

ksmooth<-function(...)
{	
	call<-match.call(stats:::ksmooth)
	ans<-stats:::ksmooth(...)
   ret<-new("ksmooth", x=ans$x, y = ans$y, call=call)

return(ret)
}


setMethod("plot", signature(x="ksmooth",y="missing"), function(x,y,...){
plot(x@x,x@y,...)
})


setMethod("lines", signature(x="ksmooth"), function(x,...){
lines(x@x,x@y,...)
})

setClass("locpoly", representation(x="numeric", y="numeric", call="call"))

locpoly<-function(...)
{
	call<-match.call(KernSmooth:::locpoly)
	ans<-KernSmooth:::locpoly(...)
	ret<-new("locpoly", x=ans$x, y=ans$y, call=call)
return(ret)
}

setMethod("plot", signature(x="locpoly",y="missing"), function(x,y,...){
plot(x@x,x@y,...)
})


setMethod("lines", signature(x="locpoly"), function(x,...){
lines(x@x,x@y,...)
})




