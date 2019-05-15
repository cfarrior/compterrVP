rm(list=ls())
setwd("c:/usr/comp_terr/compterrVP/")
source("compterr_barephys.r")

#setting up file locations
mainfile = "c:/usr/comp_terr/compterrVP/"

outfolder = paste(mainfile,"RoutputBare/",sep="")
dir.create(outfolder)
figfolder = paste(mainfile,"FiguresBare/",sep="")
dir.create(figfolder)

#better colors than standard
library(RColorBrewer)

#labels to be used in the plotting
x.lab = expression(paste(italic(x),", proportion of roots sent away"))
n.lab = expression(paste(italic(n),", number of individuals"))
fec.lab = expression(paste(italic(f),", fecundity (ind/yr)"))
r.lab = expression(paste(italic(r),", fine root biomass (g)"))
cx.lab = expression(atop(paste(italic(c),scriptstyle(x),", cost of sending roots"),"away (g/g)"))
nd.lab = expression(paste(italic(n),scriptstyle(d),", individual density (ind/d ",m^2,")"))
mort.lab = expression(paste(italic(symbol(mu)),", mortality rate (1/yr)"))
rd.lab = "fine root density, (g/dm2)"
fd.lab = "plot fecundity (ind/dm2)"
cxnot.lab = expression(paste(italic(c),scriptstyle(x),scriptstyle(0),", cost of sending roots away (g/g/dm)"))
pcompo.lab = "competitive over-investment"
Rlab = expression(paste(italic(R),", resource input rate (g",dm^2,"/yr)"))

#some plotting defaults to be used
pcx = 1
cxlab = 1.1
maxnd = 1
cols = brewer.pal(4,"Dark2")

#letterV for automatic labelling
letterV = c("A","B","C","D","E","F","G","H","I")

rundata = FALSE; #When false, suppresses running all of the anlyses but reads them from file. Fast way to edit plotting without generating data each time. 
savefig = FALSE #when true save plots to pdf in figfolder.

mort = 1

#Figure 1: cartoons to illustrate the models, produced in powerpoint.
	
{#FIGURE 2. Model 1: influence of number of individuals competing and cost of competition 

	#cx = 0, exploration of n variation
	filestem = "nvar_cx0" 
	R = 2; delta = .02; cx = 0
	comp = "ncomp"
	if(rundata)	for(n in seq(1,8)) getESS_main(comp=comp,n=n,R=R,delta=delta,cx=cx,nd=nd,init=n==1,filestem=filestem,popdynEq=FALSE)

	datac0 = read.table(paste(outfolder,filestem,".txt",sep=""),sep="\t",header=TRUE)
	datac0[datac0$fESS<0,]$fESS = 0 
	datac0[1,2] = 0 #when n = 1, x is meaningless, should be 0
		
	#cx = 1.5 (same as cxnot = 1, nd = 1), exploration of n variation
	R = 2; delta = 0.02; cr = 1.7; cx = 1.5
	filestem = "ncxvar"; plotting = FALSE; i = 0 
	if(rundata) for(n in seq(1,8,by=1)){
		i = i + 1 
		getESS_main(comp="ncomp",filestem=filestem,n=n,R=R,delta=delta,cx=cx,cxnot=cx,init=i==1,popdynEq=FALSE)
	}
	datac1 = read.table(paste(outfolder,filestem,".txt",sep=""),sep="\t",header=TRUE)
	datac1[datac1$fESS<0,]$fESS = 0 
	datac1$pcompo = (datac1$foptim-datac1$fESS)/datac1$foptim

	xplots = 2; yplots = 2; 
	if(savefig) pdf(paste(figfolder,"nvar",".pdf",sep=""),width=xplots*2.5,height=yplots*2.5)
	if(!savefig) X11(width=xplots*2.5,height=yplots*2.5)
	par(mfrow=c(yplots,xplots),cex.lab=cxlab)
	
	plot(datac0$n,datac0$r,ylab="",xlab="",ylim=c(0,1.5),pch=19,col="gray",cex=pcx)
	lines(seq(2,10),seq(2,10)*0 + (R/cr-delta/uS),col="gray")
	mtext("A",adj=0,line=2)
	mtext(r.lab,side=2,line=2.5,cex=.8)
	mtext(n.lab,side=1,line=2.5,cex=.8)
	points(datac1$n,datac1$r,pch=19,col="black")

	plot(datac0$n,datac0$x,ylab="",xlab="",ylim=c(0,1),pch=19,col="gray",cex=pcx)
	lines(seq(2,10,by=0.01),1-1/seq(2,10,by=0.01),col="gray")
	mtext("B",adj=0,line=2)
	mtext(x.lab,side=2,line=2.5,cex=.8)
	mtext(n.lab,side=1,line=2.5,cex=.8)
	points(datac1$n,datac1$x,pch=19,col="black")
	
	plot(datac0$n,datac0$fESS,ylab="",xlab="",pch=19,col="gray",cex=pcx,ylim=c(0,max(datac0$fESS)+.5))
	lines(seq(2,10),seq(2,10)*0,col="gray")	
	mtext("C",adj=0,line=2)
	mtext(fec.lab,side=2,line=2.5,cex=.8)
	mtext(n.lab,side=1,line=2.5,cex=.8)
	points(datac1$n,datac1$fESS,pch=19,col="black")
	legend(4,2.2,c("cx=0","cx=1.5","mean-field"),pch=c(19,19,NA),col=c("gray","black","gray"),lty=c(NA,NA,1),bty="n",cex=.75)
	
	plot(datac0$n,(datac0$foptim-datac0$fESS)/datac0$foptim,ylim=c(0,1.05),ylab="",xlab="",pch=19,col="gray",cex=pcx)
	lines(seq(2,10),seq(2,10)*0+1,col="gray")	
	mtext("D",adj=0,line=2)
	mtext("degree of overproliferation",side=2,line=2.5,cex=.8)
	mtext(n.lab,side=1,line=2.5,cex=.8)
	points(datac1$n,datac1$pcompo,pch=19,col="black")
	text(5,0.2,expression(paste("(",italic(f),scriptstyle(max),"-",italic(f),scriptstyle(ESS),")"/italic(f),scriptstyle(max),sep="")))
	
	if(savefig) dev.off()		
}# end Figure 2


{#Supplement figure B1. Model 1: effect of cx, given n = 3

	#influence of cx 	
	R = 2; delta = .02; n = 3; cr = 1.7
	filestem = "cxvar"; plotting = FALSE; i = 0 
	if(rundata) for(cx in seq(0,8,by=0.25)){
		i = i + 1 
		getESS_main(comp="ncomp",filestem=filestem,n=n,R=R,delta=delta,cx=cx,cxnot=cx,init=i==1,popdynEq=FALSE)
	}
	cxdata = read.table(paste(outfolder,filestem,".txt",sep=""),sep="\t",header=TRUE)
	cxdata$pcompo = (cxdata$foptim-cxdata$fESS)/cxdata$foptim

	xplots=2; yplots=2; 
	if(savefig) pdf(paste(figfolder,filestem,".pdf",sep=""),width=xplots*2.5,height=yplots*2.5)
	if(!savefig) X11(width=xplots*2.5,height=yplots*2.5)
	par(mfrow=c(yplots,xplots),cex.lab=cxlab)
	for(column in c(2,1,3,13)){ #(2) x, (1) r, (3) fESS, (13) fraction of reproduction lost due to belowground competition
		if(column==2){ y.lab=  x.lab; letter = "A";ylim=c(0,.8)}
		if(column==1){ y.lab = r.lab; letter = "B";ylim=c(0,2)}
		if(column==3){ y.lab = fec.lab; letter = "C";ylim=c(0,2)}
		if(column==13){ y.lab = pcompo.lab; letter = "D";ylim=c(0,1)}
		plot(cxdata$cx,cxdata[,column],type="n",ylab="",xlab="",ylim=ylim)
		lines(cxdata$cx,cxdata[,column],pch=19,col="black",lwd=2)
		mtext(y.lab,side=2,line=2.5,cex=.8)
		mtext("Cost of sending roots \naway (cx)",side=1,line=3,cex=.8)
		mtext(letter,adj=0,line=2)

		if(column==13) text(4.75,0.9,expression(paste("(",italic(f),scriptstyle(max),"-",italic(f),scriptstyle(ESS),")"/italic(f),scriptstyle(max),sep="")))
	}
	if(savefig) dev.off()
	
}#end Supplement figure
	

{#Figure 3. Model 2. Distance/number of competitors pairwise invasion plot. 

	#this seems to be the result no matter if you do popDynEq, r, x at ESS's or different environments, etc... (more detail in the text)
	R = 2; delta = 0.02
	aa = NULL
	r = 1; x = 0.5; cx = 1.5
	for(nrings in c(0,1,2)){
		paramsV = list(comp="hex",n=NaN,R=R,cx=cx,cxnot=NaN,cpos=cpos,delta=delta,nd=1,nrings=nrings)
		for(nrings.i in c(0,1,2)){
			inv = get_inv(r,r,x,x,nrings,nrings.i,paramsV)
			aa = rbind(aa,data.frame(nrings,nrings.i,inv))
		}
	}		

	xplots=1; yplots=1
	if(savefig) pdf(paste(figfolder,"hexdistinv",".pdf",sep=""),width=xplots*3.5,height=yplots*3.5)
	if(!savefig) X11(width=xplots*3.5,height=yplots*3.5)
	par(mfrow=c(yplots,xplots),cex.lab=cxlab)
	
	aa$inv = round(aa$inv,4)

	plot(aa$nrings,aa$nrings.i,type="n",xlab="Resident's #rings reached",ylab="Mutant's #rings reached",ylim=c(-.2,2.2),xlim=c(-.2,2.2),xaxt="n",yaxt="n")
	axis(1,c(0,1,2),c("0","1","2"))
	axis(2,c(0,1,2),c("0","1","2"))
	winners = aa[aa$inv>0,]
	points(winners$nrings,winners$nrings.i,pch = "+",cex=2)
	losers = aa[aa$inv<0,]
	points(losers$nrings, losers$nrings.i, pch="-",cex=2)
	ok = aa[aa$inv==0,]
	points(ok$nrings,ok$nrings.i,pch="O")
	abline(0,1,lty=2)	
	if(savefig) dev.off()

}# end Figure 3

{#Figure 4, Figure C2,C3,C4. Model 2. Predictions across environmental gradients, under various assumptions of the cost of moving roots through soil (cxnot) and leakage (delta). 

	#fine = TRUE #use fine=FALSE when wanting to run faster and check coarse effects. 
	
	lwd = 1; reps = 3; 
	for(deltaD in c(0,0.02)){ #default values of delta
		for(cxnotD in c(0,1.2)){ #default values of cxnot
			fileadd = ""
			if(fine) fileadd = "FINE"
			comp = "hex"; ofile = paste(outfolder,"hexcomp",fileadd,dim(cpos)[1],"_",cxnotD,"_",deltaD*1000,"_",sep="")
			
			if(!fine) nums = 4; if(fine) nums = 20
			Rmin= 0.1; Rmax = 2.25
			RV = seq(Rmin,Rmax,by=(Rmax-Rmin)/(nums-1))
				
			if(rundata){
				a = NULL; delta = deltaD; cxnot = cxnotD; ii = 0 
				#Do not run RV if both cxnot and delta are 0 because in this case there is no equilibrium of individual density. The full tragedy of the commons in this case sends plants to a point of no reproduction.
				if(cxnotD+deltaD!=0) for(R in RV){
					out = try(runhexdistcomp(R,delta,cxnot,comp=comp))
					if(is.numeric(out[[1]])){ 
						ii = ii + 1
						go = data.frame(R,delta,cxnot,out)
						print(go)
						write.table(go,paste(ofile,"Rvar.txt",sep=""),sep="\t",row.names=FALSE,col.names=ii==1,append=ii!=1)
						a = rbind(a,go)
					}
				}
				if(cxnotD+deltaD!=0) Rvar = a
				
				R = 2; delta=deltaD; a = NULL; ii = 0 
				
				cxmin = 0; cxmax = 3.5
				cxnotV = seq(cxmin,cxmax,by=(cxmax-cxmin)/(nums-1))
				for(cxnot in cxnotV){
					out = try(runhexdistcomp(R,delta,cxnot,comp=comp))
					if(is.numeric(out[[1]])){
						ii = ii + 1
						go = data.frame(R,delta,cxnot,out)
						print(go)
						write.table(go,paste(ofile,"cxvar.txt",sep=""),sep="\t",row.names=FALSE,col.names=ii==1,append=ii!=1)
						a = rbind(a,go)
					}
				}
				cxvar = a
				
				R = 2; cxnot = cxnotD; a = NULL; ii = 0 
				dmin = 0; dmax = 0.4
				deltaV = seq(dmin,dmax,by=(dmax-dmin)/(nums-1))
				for(delta in deltaV){
					out = try(runhexdistcomp(R,delta,cxnot,comp=comp))
					if(is.numeric(out[[1]])){ 
						ii = ii + 1
						go = data.frame(R,delta,cxnot,out)
						print(go)
						write.table(go,paste(ofile,"delvar.txt",sep=""),sep="\t",row.names=FALSE,col.names=ii==1,append=ii!=1)
						a = rbind(a,go)	
					}
				}
				delvar = a
			}
			
			if(deltaD+cxnotD!=0) Rvar = read.table(paste(ofile,"Rvar.txt",sep=""),header=TRUE,sep="\t")
			if(deltaD+cxnotD==0) Rvar = data.frame(R=RV,r=0,x=0,nd=0,f=0,foptim=1)
			cxvar = read.table(paste(ofile,"cxvar.txt",sep=""),header=TRUE,sep="\t")
			delvar = read.table(paste(ofile,"delvar.txt",sep=""),header=TRUE,sep="\t")

			#this version of the model, without investment in height is only valid for nd<1. 
			Rvar = Rvar[Rvar$nd<1,]; Rvar = Rvar[order(Rvar$R),]
			cxvar = cxvar[cxvar$nd<1,]; cxvar = cxvar[order(cxvar$cxnot),]
			delvar = delvar[delvar$nd<1,]; delvar = delvar[order(delvar$delta),]
			
			xplots = 3; yplots = 1
			if(!savefig) X11(width=xplots*2.5,height=yplots*2.5)
			if(savefig) pdf(paste(figfolder,comp,"_cxnot",cxnotD*1000,"_delta",deltaD*1000,".pdf",sep=""),width=xplots*2.5,height=yplots*2.5)
			par(mfrow=c(yplots,xplots),cex.lab=1.1,oma=c(0,0,0,2),mar=c(4,4,4,2)+0.1)
			
			for(var in c("R","cx","del")){
				if(var=="R") {data = Rvar; x =data$R; xlab=expression(paste(italic(R),", resourcee input (g ",dm^-2,yr^-1,")"));letter="A"}
				if(var=="cx"){data = cxvar; x =data$cxnot; xlab=expression(paste(italic(c),scriptstyle(x),scriptstyle(0),", (",dm^-1, yr^-1,")"));letter="B"}
				if(var=="del"){data = delvar; x =data$delta; xlab=expression(paste(symbol(Delta),", leakage loss (",yr^-1,")"));letter="C"}
					
				plot(x,data$x,ylim=c(0,1.1),type="n",yaxt="n",xaxt="n",xlab="",ylab="")
				if(var=="del") mtext(expression(paste(italic(r)," (g)")),side=4,line=2.5,cex=.75)
				lines(x,data$x,col="red",lwd=1,lty=3)
				lines(x,data$nd,lwd=lwd,lty=1)
				axis(side=2,at=c(0,.5,1),cex=.5)
				x1=0.5;y1=1.1
				if(var=="R") legend(x1,y1,legend=c(expression(paste(italic(r),scriptstyle(ESS),sep="")),expression(paste(italic(x),scriptstyle(ESS),sep="")),expression(paste(italic(n),scriptstyle(d),sep=""))),
				bty="n",lwd=c(1,1),col=c("darkgreen","red","black"),lty=c(2,3,1),cex=.75)				

				par(new=TRUE)
			
				ymax = 4.5
				if(cxnotD ==0) ymax = 7
				plot(x,data$r,col="brown",pch=19,ylim=c(0,ymax),xlab="",type="n",ylab="",yaxt="n",xaxt="n")
				mtext(xlab,side=1,cex=.75,line=3)
				if(var=="R") mtext(expression(paste(italic(x)," and ",italic(n),scriptstyle(d))),side=2,line=2.5,cex=.75)

				axis(side=1,cex=.5)
				lines(x,data$r,col="darkgreen",lwd=lwd,lty=2)
				axis(side=4,cex=.5)	
				mtext(letter,adj=0,line=2)				
			}	
			if(savefig) dev.off()
		}
	}
}#end Figure 4, Figure C2,C3,C4


{#Figure C5. Comparison of the ESS and population-level optimum of allocation to roots. Qualitative patterns are the same for r and r*nd as plotted. 
	deltaD = 0.02; cxnotD = 1.2
	fileadd = ""
	if(fine) fileadd = "FINE"
	comp = "hex"; ofile = paste(outfolder,"hexcomp",fileadd,dim(cpos)[1],"_",cxnotD,"_",deltaD*1000,"_",sep="")
	
	Rvar = read.table(paste(ofile,"Rvar.txt",sep=""),header=TRUE,sep="\t")
	cxvar = read.table(paste(ofile,"cxvar.txt",sep=""),header=TRUE,sep="\t")
	delvar = read.table(paste(ofile,"delvar.txt",sep=""),header=TRUE,sep="\t")
	Rvar = Rvar[Rvar$nd<1,]; Rvar = Rvar[order(Rvar$R),]
	cxvar = cxvar[cxvar$nd<1,]; cxvar = cxvar[order(cxvar$cxnot),]
	delvar = delvar[delvar$nd<1,]; delvar = delvar[order(delvar$delta),]
		
	xplots = 3; yplots = 1
	if(!savefig) X11(width=xplots*2.5,height=yplots*2.5)
	if(savefig) pdf(paste(figfolder,"rcompoptim",cxnotD*10,deltaD*1000,".pdf",sep=""),width=xplots*2.5,height=yplots*2.5)
	par(mfrow=c(yplots,xplots),cex.lab=1.1,oma=c(0,0,0,2),mar=c(4,4,4,2)+0.1)
	
	for(var in c("R","cx","del")){
		if(var=="R"){data = Rvar; x =data$R; xlab=expression(paste(italic(R),", rainfall (dm ",yr^-1,")"));letter="A"}
		if(var=="cx"){data = cxvar; x =data$cxnot; xlab=expression(paste(italic(c),scriptstyle(x),scriptstyle(0),", root extension cost"));letter="B"}
		if(var=="del"){data = delvar; x =data$delta; xlab=expression(paste(symbol(Delta),", leakage loss (",yr^-1,")"));letter="C"}
		
		plot(x,data$r*data$nd,col="brown",pch=19,ylim=c(0,1.5),xlab="",type="n",ylab="",yaxt="n",xaxt="n")
		mtext(xlab,side=1,cex=.75,line=3)
		lines(x,data$r*data$nd,col=cols[2])
		lines(x,data$roptim*data$nd,col=cols[3])
		if(var=="R"||var=="cx") mtext(expression(paste("fine root biomass (g)")),side=2,line=2)
		axis(side=1,cex=.5)
		axis(side=2,cex=.5)
		
		if(var=="R") legend(0.5,1.5,legend=c(expression(paste(italic(r),scriptstyle(ESS),"*",italic(n),scriptstyle(d),sep="")),expression(paste(italic(r),scriptstyle(OPTIM),"*",italic(n),scriptstyle(d),sep=""))),bty="n",col=c(cols[2],cols[3]),lwd=c(1,1),lty=c(1,1),cex=.75)
		
	}	
	if(savefig) dev.off()		
}

{#Figures B1, B2. Model 1: ESS and CSS verification plots. 

	#Model 1: appendix graphs of 2d optimization
	savefig = TRUE
	filestem = "nvar_cx0" 
	datac0 = read.table(paste(outfolder,filestem,".txt",sep=""),sep="\t",header=TRUE)
	ESSCSSverificationplots(datac0,comp="ncomp",var="nvarc0")	
	
	filestem = "ncxvar"
	datac1 = read.table(paste(outfolder,filestem,".txt",sep=""),sep="\t",header=TRUE)
	ESSCSSverificationplots(datac1,comp="ncomp",var="nvarc2")	

}#end Figure B1	


{#Figure B3, B4, B5, B6; ESS and CSS and population size stability verification for Figures 4, plus figures for C2,C3,C4 (C appendix figure checks produced here, but not printed in the supplement)
	fine = TRUE 
	
	lwd = 1; reps = 3; 
	for(deltaD in c(0,0.02)){
		for(cxnotD in c(0,1.2)){
			fileadd = ""
			if(fine) fileadd = "FINE"
			comp = "hex"; ofile = paste(outfolder,"hexcomp",fileadd,dim(cpos)[1],"_",cxnotD,"_",deltaD*1000,"_",sep="")
			
			if(!fine) nums = 4; if(fine) nums = 20
			Rmin= 0.5; Rmax = 2.25
			RV = seq(Rmin,Rmax,by=(Rmax-Rmin)/(nums-1))

			if(deltaD+cxnotD!=0) Rvar = read.table(paste(ofile,"Rvar.txt",sep=""),header=TRUE,sep="\t")
			if(deltaD+cxnotD==0) Rvar = data.frame(R=RV,r=0,x=0,nd=1,delta=deltaD,cx=cxnotD,cxnot=cxnotD, nrings.i =1, nrings=1)
			cxvar = read.table(paste(ofile,"cxvar.txt",sep=""),header=TRUE,sep="\t")
			cxvar = cxvar[!is.na(cxvar$inv),]
			delvar = read.table(paste(ofile,"delvar.txt",sep=""),header=TRUE,sep="\t")
			
			presplot = TRUE; 
			rV = seq(0.1,10,by=0.1)
			for(var in c("R","cx","del")){
				if(var=="R") data = Rvar
				if(var=="cx") data = cxvar
				if(var=="del") data = delvar

				data$n  = 7; data$nrings = 1
				
				if(length(unique(data$nrings.i))>1) data = data[data$nrings.i==1,] #hackey way of keeping it from repeating when x=0 is the answer (due to how the code returns the distance ess when xESS=0).
				
				data = data[data$nd<1,]
				
				#optional step for speed, subsetting the data to run ESS/CSS checks
				if(dim(data)[1]>4){
					data = rbind(data[1,],data[round(dim(data)[1]/3),],data[round(dim(data)[1]/3*2),],data[dim(data)[1]-1,])	
				}

				if(dim(data)[1]>0) ESSCSSverificationplots(data,comp=comp,var=var,fileadd=paste(fileadd,cxnotD*10,"_",deltaD*100,sep=""))	
			}
		}						
	}
}#end Figure B3,B4,B5,B6
