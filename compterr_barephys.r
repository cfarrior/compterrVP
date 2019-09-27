#C. Farrior - last updated Sept 26, 2019 
#cfarrior@gmail.com

#This file holds the meat of the biology in the model.
#Fitness functions, photosynthesis functions, etc..

#needed libraries and other files
library(plotrix) #used for plotting vectorfields
source("ess2d.r")
source("hexbuild.r")

#defaults for debugging, etc... 
plotting = TRUE #watch all of the optimization happen. This makes way too many plots when not debugging. 
filestem = "test" #default output filestem
init = TRUE #write to a new file 
project = "compterr"

#vectors giving the range and accuracy of the optimized traits
xV =  seq(0,1,by=0.01)
rV =  seq(0.1,20,by=0.1)

#competition parameters and costs
comp = "ncomp" #competition/model type "ncomp" or "hex"
n = 7 #number of competitors (including the target individual)
nd = 1 #individual density 
cxnot = 2 #cost per unit root per unit distance to nearest neighbors of all roots sent outside of "home"
uS = 1	#uptake efficiency of fine roots (resource/root)
cf = 1 #cost of producing a new individual
cr = 1.7 #cost of roots, independent of placement

################
get_cx = function(cxnot,nd,n){
################
#determine the distance-related cost of sending roots, given a standard cxnot and density of individuals. 
	if(n==7){
		hexa = (1/nd*2/(3*3^.5))^.5 #the length of a side of hexagon if they are equally sized and connected perfectly given nd density of individuals. 
		bdist = hexa*3^.5 #the distance between two closest hexagon midpoints
		cx = cxnot*bdist^2 		
		return(cx)
	}
	if(n==2) return(cxnot*n/pi/nd)
}#end get_cx


#some more defaults
cx = get_cx(cxnot,nd,n)


nrings = NaN

reps = 7
presplot = FALSE
gohex = TRUE

if(gohex){
	posV = get_hex(4); pos = posV[[1]]; neighborsM = posV[[2]]
	#pos columns: (1) x (2) y (3) unum (4) inv/res
	#where x,y are coordinates in eucliddean space; unum = unique individual number; inv/res 
	cpos = cbind(pos,0)
	#put the invader in the middle of the invididuals 
	inv = closest(mean(cpos[,1]),mean(cpos[,2]));  cpos[cpos[,3]==inv[,3],4]=1
}

#grid used for plotting of competeem2d for ESS, CSS verification plots
var1Vgo = seq(0.1,5,by=.5); var2Vgo = seq(0,1,by=0.1)
rxcomp_Mtemp = NULL; var1tick = 0; 
for(var1 in var1Vgo){
	var1tick = var1tick+1; var2tick = 0 
	for(var2 in var2Vgo){
		var2tick = var2tick + 1
		rxcomp_Mtemp  = rbind(rxcomp_Mtemp ,c(var1,var2,var1tick,var2tick))
	}
}

#needed for plotting. Let me know if there is a better way! 
get_os <- function() {
  if (.Platform$OS.type == "windows") { 
    "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "mac" 
  } else if (.Platform$OS.type == "unix") { 
    "unix"
  } else {
    stop("Unknown OS")
  }
}
os = get_os()

##############################
getParamsV = function(){
#organizing function to collect parameter values in global environment. Only used rarerly in code, but useful for debugging. 
##############################

	paramsV = list(comp=comp,n=n,R=R,cx=cx,cxnot=cxnot,cpos=cpos,delta=delta,nd=nd,nrings=nrings,mort=mort)
	
	if(comp=="hex"){
		cx = get_cx(cxnot,nd,n)
		paramsV$cx = cx; paramsV$n = n; nd = paramsV$nd
	}
		
	return(paramsV)
}#end getParamsV

##############################
getESS_main = function(comp,n=NaN,R=1,delta=0,cx=NaN,nd=1,cpos=NaN,init=NaN,filestem="test",cxnot=NaN,popdynEq=FALSE,mort=NaN,write=TRUE,nrings=NaN,justParamsV=FALSE){
#comp: competition type, "ncomp" is model 1, n individual competition model and "hex" is honeycomb competition
#n: used in "ncomp" only, the number of individuals competing. 
#R: the constant resource input rate
#cx: cost per unit fine root of sending that unit of fine root to a neighbor, the n-1 neighbors in "ncomp" and the 6 closest neighbors in "hex"
#nd: individual density, equals 1 for "ncomp"
#cpos: (used only if comp=="hex") matrix with the x,y positions of each of the individuals in competition columns: (1) x (2) y (3) unique individual number (4) invader (1) or resident (0) type
#init: logical, if TRUE overwrites file, if FALSE appends to file for saving output; if NaN turns into TRUE if no file exists yet FALSE if one does
#popdynEq: if TRUE find the nd that equilibrates (by changing cx) 
#nrings: #rings of neighbors to compete over, only matters if in comp = "hex"
##############################
	rinit = 1; xinit = 0.5
	
	#list of parameters needed to be passed through all functions
	paramsV = list(comp=comp,n=n,R=R,cx=cx,cxnot=cxnot,cpos=cpos,delta=delta,nd=nd,nrings=nrings,mort=mort)
	
	if(is.na(cx)){ 
		cx = get_cx(cxnot,nd,n)
		paramsV$cx = cx; paramsV$n = n; nd = paramsV$nd
	}
	if(justParamsV) return(paramsV)
	
	if(popdynEq){	
		rpdINIT=NaN; xpdINIT = NaN
		if(!is.na(nrings)) if(nrings==0){
			rpdINIT = PIP_ESS_telescope(var="var1",rV,1,0,paramsV)
			xpdINIT = 0 
		}
		if(is.na(xpdINIT)){
			essV = getESS_2D(rinit,xinit,rV,xV,reps,paramsV)		
			rpdINIT = essV[[1]]; xpdINIT = essV[[2]]
		}

		fV = get_f(rpdINIT,xpdINIT,paramsV)
		f = fV$f;

		undermaxnd = TRUE
		if(f>=mort){
			nd = 1 #nd>=1, not considered in this model. See BioRxiv model for an implementation of height competition investment when the canopy closes. 
			fV = get_f(rpdINIT,xpdINIT,paramsV)
			result = list(cx=paramsV$cx,nd=nd,r=rpdINIT,x=xpdINIT,f=fV$f,write=FALSE)
			undermaxnd = FALSE
		}

		if(rpdINIT>=max(rV)-.5) print ("error")
				
		if(undermaxnd){
			nd = try(uniroot(zeroPopDyn,interval=c(0.001,maxnd*2),paramsV=paramsV)[[1]],silent=TRUE)
			if(!is.numeric(nd)) nd = try(uniroot(zeroPopDyn,interval=c(0.05,maxnd*2),paramsV=paramsV)[[1]],silent=TRUE)
			if(!is.numeric(nd)) nd = try(uniroot(zeroPopDyn,interval=c(0.1,maxnd*2),paramsV=paramsV)[[1]],silent=TRUE)
			if(!is.numeric(nd)){ nd = 0; print("nd find error")}
			paramsV$nd = nd 
			paramsV$cx = get_cx(cxnot,nd,n)
			
			#check that it's negative and not positive density dependence.
			bb = NULL
			ndV = seq(max(0.1,nd-.1),nd+.1,by=0.025)
			if(min(ndV)>=nd) ndV = seq(nd-.2*nd,nd+.2*nd,by=.05*nd)
			for(ndgo in ndV){
				fminusm = zeroPopDyn(ndgo,paramsV) #within zeroPopDyn the n is updated with ndgo
				bb = rbind(bb,c(ndgo,fminusm))
			}
			if(os!="mac") X11()
			if(os=="mac") quartz()
			plot(bb[,1],bb[,2],main="",xlab="nd",ylab="f-m")
			lines(bb[,1],bb[,2])
			abline(0,0)
			lines(c(nd,nd),c(-100,100))
			
			r=NaN; x = NaN
			if(!is.na(nrings)) if(nrings==0){
				r = PIP_ESS_telescope(var="var1",rV,mean(rV),0,paramsV)
				x = 0 
			}
			if(!is.na(cxnot)) if(cxnot==0){
				r = PIP_ESS_telescope(var="var1",rV,mean(rV),1-1/n,paramsV)
				x = 1-1/n 
			}
			if(is.na(x)){
				essV = getESS_2D(rinit,xinit,rV,xV,reps,paramsV)		
				r = essV[[1]]; x = essV[[2]]
			}
		}
	}
	
	if(popdynEq) if(!undermaxnd){r=rpdINIT;x=xpdINIT}
	
	if(!popdynEq){
		essV = getESS_2D(rinit,xinit,rV,xV,reps,paramsV)		
		r = essV[[1]]; x = essV[[2]]
	}
	if(r<=min(rV)) x = 0 

	fV = get_f(r,x,paramsV)
	f = fV$f.i 

	roptim = NaN; foptim=NaN
	roptim = try(optim(par=1,fn=get_f_hex,paramsV=paramsV,roptim=TRUE,control=list(fnscale=-1),lower=0),silent=TRUE)
	if(!is.numeric(roptim[[1]])) roptim = optim(par=1,fn=get_f_n,paramsV=paramsV,roptim=TRUE,control=list(fnscale=-1),lower=0)

	foptim = roptim$value	
	roptim = roptim$par

	if(plotting){
		competeem2d(rxcomp_Mtemp,.1,.01,paramsV,plotting=TRUE)
		points(r,x,col="red",pch=19)
	}
	
	if(write){ 
		if(is.na(init)){
			filetry = try(read.table(paste(outfolder,filestem,".txt",sep=""),header=TRUE))
			if(is.numeric(filetry)[[1]]) init=FALSE
			if(!is.numeric(filetry)[[1]]) init=TRUE
		}
		write.table(data.frame(r,x,fESS=f,comp=paramsV$comp,n=paramsV$n,cx=paramsV$cx,R=paramsV$R,delta=paramsV$delta,nd=paramsV$nd,foptim=foptim,roptim=roptim,nrings=nrings),paste(outfolder,filestem,".txt",sep=""),append=!init,sep="\t",row.names=FALSE,col.names=init)
	}	
		
	return(list(r=r,x=x,f=f,cx=paramsV$cx,nd=paramsV$nd,roptim=roptim,foptim=foptim))
}#end getESS_main 

##############################
zeroPopDyn = function(nd,paramsV,r=NaN,x=NaN){	
#returns zero if population dynamics are stable
##############################
	rinit = 1; xinit = 0.5
	
	paramsV$nd = nd 
	
	paramsV$cx = get_cx(paramsV$cxnot,nd,paramsV$n)
	
	if(is.na(r)) if(is.na(x)){
		if(!is.na(paramsV$nrings)) if(paramsV$nrings==0){
			r = PIP_ESS_telescope(var="var1",rV,1,0,paramsV)
			x = 0 
		}
	}	
	if(is.na(x)){
		essV = getESS_2D(rinit,xinit,rV,xV,reps,paramsV)		
		r = essV[[1]]; x = essV[[2]]
	}

	f = get_f(r,x,paramsV)$f
	
	return(f-paramsV$mort)
}#end zeroPopDyn

##############################
ESSCSSverificationplots = function(data,comp="hex",var,fileadd=""){
#function used to check that ESSmain returns both a result that is both an ESS and a CSS
##############################
	
	for(i in seq(1,dim(data)[1])){
		line = data[i,]
		R = line$R; delta=line$delta; nrings = line$nrings; nd = line$nd
		if(comp =="hex") cxnot = line$cxnot
		if(comp == "ncomp"){
			cx = line$cx
			n = line$n
			nrings=1
		}
		paramsV = list(comp=comp,n=n,R=R,cx=cx,cxnot=cxnot,cpos=cpos,delta=delta,nd=nd,nrings=nrings,mort=mort)
		if(comp=="hex"){
			cx = get_cx(cxnot,nd,n)
			paramsV$cx = cx; paramsV$n = n; nd = paramsV$nd
		}
		if(comp=="hex") print(paramsV$cx)
		if(comp=="ncomp") print(paramsV$n)


		if(!savefig) 
		if(!savefig){
			if(os!="mac") X11(width=8,height=3)
			if(os=="mac") quartz(width=8,height=3)
		}
			
		if(savefig) pdf(paste(figfolder,"ESSvar",var,fileadd,i,".pdf",sep=""),width=8,height=3)
		par(mfrow=c(1,3))

		#plots to show ESS (left panels)
		test_if_ess(line$r,line$x,rV,xV,paramsV,main="")
		if(comp=="hex") mtext(paste("Model 2: R=",round(line$R,1),", ",
			"cx0=",round(cxnot,2),", ",
			"Delta=",round(line$delta,3),sep=""),adj=0,line=2)
		if(comp=="ncomp") mtext(paste("Model 1: ",
			"cx=",round(cx,2),", n=",n,sep=""),adj=0,line=2)

		#plot to show convergence stability
		competeem2d(rxcomp_Mtemp,.1,.01,paramsV,plotting=TRUE,presplot=TRUE)
		points(line$r,line$x,col="green",pch=19)

		#plot to show population dynamic equilibrium
		if(comp=="hex"){
			bb = NULL
			ndV = seq(max(0.1,nd-.1),nd+.1,by=0.025)
			if(min(ndV)>=nd) ndV = seq(nd-.2*nd,nd+.2*nd,by=.05*nd)
			for(ndgo in ndV){
				plotting=FALSE
				fminusm = zeroPopDyn(ndgo,paramsV) #within zeroPopDyn the n is updated with ndgo
				bb = rbind(bb,c(ndgo,fminusm))
			}
			plot(bb[,1],bb[,2],xlab=expression(paste(italic(n),scriptstyle(d),", Site individual density")),ylab=expression(paste(italic(f-m),", Population growth rate")))
			lines(bb[,1],bb[,2])
			abline(0,0,lty=2)
			lines(c(nd,nd),c(-100,100),col="green")
		}
		if(savefig) dev.off()
	}		
}#end ESSCSSverificationplots


##############################
runhexdistcomp = function(R,delta,cxnot,filestem="hextesting",comp="hex"){
#Function used to find the competitive dominant allocation of fine roots through space (how many rings of competitors to engage with)..
##############################

	#This sets the distance of the resident's fine root spread to the first ring of competitors (neighbors will be found in "hex" space <1.1 hex units away)
	nrings.r = 1
	
	#This is not used in the paper, but the ncomp model with n = 7 and population dynamics is a good approximation of the much slower "hex" model. 
	if(comp=="ncomp"){
		essV = getESS_main(comp="ncomp",n=7,R=R,delta=delta,cx=NaN,nd=1,cpos=cpos,init=NaN,filestem=filestem,cxnot=cxnot,popdynEq=TRUE,mort=mort,write=TRUE,nrings=nrings.r) 
		return(essV)
	}

	if(comp=="hex"){
		essV = getESS_main(comp="hex",n=7,R=R,delta=delta,cx=NaN,nd=1,cpos=cpos,init=NaN,filestem=filestem,cxnot=cxnot,popdynEq=TRUE,mort=mort,write=TRUE,nrings=nrings.r) 
		
		aa = NULL
		r = essV[[1]]; x = essV[[2]]; cx = essV$cx; nd = essV$nd; 
		for(nrings in c(nrings.r)){
			paramsV = list(comp="hex",n=7,R=R,cx=cx,cxnot=NaN,cpos=cpos,delta=delta,nd=nd,nrings=nrings)
			
			for(nrings.i in c(0.1,1.1,2.5)){
				inv = get_inv(r,r,x,x,nrings,nrings.i,paramsV)
				aa = rbind(aa,data.frame(nrings,nrings.i,inv,essV))
			}
		}
		#this returns the best strategies in competition with the resident that competes with individuals in the closest ring. 
		return(aa[aa[,3]==max(aa[,3]),])
		#note: if staying at home is optimal xESS will equal zero and nrings.r =1 will not matter. 
	}
}#end runhexcompdist
	

##############################
get_f = function(r,x,paramsV){
#function that only shuttles the fitness call to the specific fitness functions of the separate models
##############################	
		
	if(paramsV$comp=="ncomp") fV = get_f_n(r.r=r,r.i=r,x.r=x,x.i=x,paramsV=paramsV)
	if(paramsV$comp=="hex") fV = get_f_hex(r.r=r,r.i=r,x.r=x,x.i=x,paramsV=paramsV)

	return(fV)
}#end get_f

##############################
get_inv = function(var1.r,var1.i,var2.r,var2.i,nrings.r=NaN,nrings.i=NaN,paramsV){
#a function called on by the ess2d routines
#takes in resident and invader values of the two traits and returns the difference between the invaders fitness and the resident fitness
##############################	

	r.r = var1.r; r.i = var1.i
	x.r = var2.r; x.i = var2.i

	if(paramsV$comp=="hex") if(is.na(nrings.r)) nrings.r = nrings.i = paramsV$nrings

	if(paramsV$comp=="ncomp") fV = get_f_n(r.r=r.r,r.i=r.i,x.r=x.r,x.i=x.i,paramsV=paramsV)
	if(paramsV$comp=="hex") fV = get_f_hex(r.r=r.r,r.i=r.i,x.r=x.r,x.i=x.i,nrings.r=nrings.r,nrings.i=nrings.i,paramsV=paramsV)
	
	return(fV$f.i - fV$f)
}#end get_inv 



##############################	
get_f_n = function(ropt=NaN,r.r=NaN,r.i=NaN,x.r=NaN,x.i=NaN,paramsV,roptim=FALSE){
#fitness function for the "ncomp", n individual competition model 
#takes in resident and invader values of traits r and x as well as all of the content in paramsV
##############################	
	if(roptim){
		r.i = r.r = ropt
		x.r = x.i = 0 
	}

	#load content of paramsV 
	n = paramsV$n; R = paramsV$R; cx = paramsV$cx; delta = paramsV$delta; nd = paramsV$nd;

	abovespace = min(1/maxnd,1/nd)
	belowspace = 1/nd
	
	#root amounts in locations for (n-1) resident individuals and 1 invading individual: 
	if(n!=1){
		r.tr = r.r*(1-x.r) + r.i*x.i/(n-1) + (n-2)/(n-1)*r.r*x.r
		r.ti = r.i*(1-x.i) + r.r*x.r
	}
	if(n==1){r.tr = r.r; r.ti=r.i}
	
	s.i = (R*belowspace)/(uS*r.ti + delta)
	s.r = (R*belowspace)/(uS*r.tr + delta)	

	if(n!=1){
		AW.i = ((1-x.i)*uS*r.i*s.i + uS*x.i*r.i*s.r)
		AW.r = ((1-x.r)*uS*r.r*s.r + (n-2)*uS*x.r/(n-1)*r.r*s.r + x.r/(n-1)*uS*r.r*s.i)
	}
	if(n==1){
		AW.i = (uS*r.i*s.i)
		AW.r = (uS*r.r*s.r)
	}
	if(r.i == 0) AW.i = 0; if(r.r ==0) AW.r = 0 
	
	f.i = 1/cf*(AW.i - cr*r.i - cx*x.i*r.i)
	f.r = 1/cf*(AW.r - cr*r.r - cx*x.r*r.r)
	
	if(roptim) return(f.i)
	return(list(f.i=f.i,f=f.r,inv=f.i>f.r,leaves=NaN))
}#end get_f_n

##################################################	
get_f_hex = function(ropt=NaN,r.r=NaN,r.i=NaN,x.r=NaN,x.i=NaN,nrings.r=NaN,nrings.i=NaN,paramsV,roptim=FALSE){
#fitness function for the "honeycomb" competition model 
#takes in resident and invader values of traits r and x as well as all of the content in paramsV
##################################################	
	if(roptim){
		r.i = r.r = ropt
		x.r = x.i = 0 
		nrings.r = nrings.i = 0 
	}
 	
	#load content of paramsV 
	R=paramsV$R; cx = paramsV$cx; nd = paramsV$nd; delta = paramsV$delta; cpos = paramsV$cpos
	
	if(is.na(nrings.r)) nrings.r=nrings.i=paramsV$nrings
	
	#nrings: this is in hex space for simplicity and not related to nd. when 1.1. neighbors are the closest six neighbors. 
	
	#calculate root density landscape
	r.M = cbind(cpos,0*seq(1,dim(cpos)[1])) #a matrix with columns (1) x (2) y (3) unum (4) inv/res (5) root density
	for(i in seq(1,dim(cpos)[1])){
		target = cpos[i,,drop=FALSE]
		
		if(target[,4]==1){ nrings = nrings.i; r = r.i; x = x.i }
		if(target[,4]==0){ nrings = nrings.r; r = r.r; x = x.r }
		
		if(nrings==0) hexdist = 0; if(nrings==1) hexdist = 1.1; if(nrings==2) hexdist = 2.5
		neighbors = get_neighbors(target,cpos,hexdist)
		
		#add the staying roots. 
		r.M[i,5] = r.M[i,5] + r*(1-x)

		#add the leaving roots to their spots. 
		Nneigh = dim(neighbors)[1]
				
		for(n.num in neighbors[,3]) r.M[n.num,5] = r.M[n.num,5] + r*x/Nneigh
	}
		
	#calculate the individual fitnesses 
	abovespace = min(1/maxnd,1/nd)
	belowspace = 1/nd

	f.M = cbind(cpos,NA*seq(1,dim(cpos)[1]))
	for(i in seq(1,dim(r.M)[1])){
		target = r.M[i,,drop=FALSE]

		if(target[,4]==0){r=r.r;x=x.r;nrings=nrings.r}
		if(target[,4]==1){r=r.i;x=x.i;nrings=nrings.i}
		
		if(nrings==0) hexdist = 0; if(nrings==1) hexdist = 1.1; if(nrings==2) hexdist = 2.5		
		neighbors = get_neighbors(target,r.M,hexdist)
		
		homeroots = r*(1-x); awayroots = r*x
				
		wtarget = uS*homeroots/(uS*target[,5]+delta)*R*belowspace
		
		if(nrings==0) hexdist = 0; if(nrings==1) hexdist = 1.1; if(nrings==2) hexdist = 2.5
		neighbors = neighbors[neighbors[,dim(neighbors)[2]]<hexdist,,drop=FALSE]

		if(dim(neighbors)[1]>0) for(j in seq(1,dim(neighbors)[1])){
			wtarget = wtarget + uS*awayroots/dim(neighbors)[1]/(uS*neighbors[j,5]+delta)*R*belowspace
		}
		
		f.M[i,5] = (1/cf)*(wtarget - cr*r - nrings^2*cx*x*r)
	}
	
	fcomp = mean(f.M[f.M[,4]==0,5]) #average fitness of all residents
	f.i = f.M[f.M[,4]==1,5] #fitness of invader
	#answers might be slightly different if you have to beat just your closest neighbors than all of the individuals on the whole torus
	
	target = f.M[f.M[,4]==1,,drop=FALSE]
			
	if(roptim) return(f.i)
	return(list(f.i=f.i,f=fcomp,inv=f.i>fcomp,leaves=NaN))

}#end get_f_hex




