#C. Farrior - last updated Dec 20, 2018 
#cfarrior@gmail.com

library(plotrix) #used for plotting vectorfields
project = NaN

#Methods for finding a 2 dimensional ESS (var1 and var2)
#This code requires a function ("get_inv") which takes in resident and invader strategies and returns whether the invader outperforms the resident.

#the necessary get_inv function is of the form: 
#get_inv(var1.r,var1.i,var2.r,var2.i,paramsV=paramsV)
#where var1 and var2 are the variables being "ESS'ed"; .r indicates the resident value; .i = indicates the invader value. 
#paramsV includes all of the parameters that you may need to execute get_inv(). this vector is passed throughout the functions in here and only openned back in your get_inv code.
#get_inv should return a positive value if the invader succeeds, 0 if it is neutral and negative if it does not succeed (usually: f.i - f.r)

##############################
getESS_2D = function(var1init,var2init,var1V,var2V,initreps=7,paramsV){
#var1init = initial guess for var1
#var2init = initial guess for var2
#var1V = vector with a min and max that encompass the range within which to look for the ess of var1, and the spacing between elements the desired accuracy of the result
#var2V = same for var2 
#initreps = the number of times to try to find the ESS by optimizing one var at a time
#paramsV = a list of parameters that are needed to compute fitness, specific to specific models
##############################
	
	if(initreps<3) initreps = 3
	
	#initialize the variables
	var1 = var1init; var2 = var2init 
	
	#start with trying to find the ESS by finding one and then the other 
	reptrack = NULL
	for(rep in seq(1,initreps)){
		var1 = PIP_ESS_telescope(var="var1",var1V,var1,var2,paramsV)
		var2 = PIP_ESS_telescope(var="var2",var2V,var1,var2,paramsV)
		reptrack = rbind(reptrack,data.frame(var1,var2))
	}
	essV = list(var1,var2)
	
	#decide whether the above process has converged on an set of strategies that make up an ESS
	#take the last 3 itterations
	last = reptrack[seq(dim(reptrack)[1]-3,dim(reptrack)[1]),]

	#if the sd is greater than the accuracy then go2d = TRUE and we will go through 2 dimensional optimizing (optimizing both var1 and var2 at once, not one at a time)
	sdvar1 = sd(last$var1); sdvar2 = sd(last$var2)
	go2d = FALSE
	if(sdvar1>var1V[2]-var1V[1]) go2d = TRUE
	if(sdvar2>var2V[2]-var2V[1]) go2d = TRUE
	
	if(go2d){
		print("going to 2d optimization")
		essV = PIP_2D_telescope(var1V,var2V,paramsV)
	}
	
	return(essV)
}


##############################
PIP_ESS_telescope = function(var,V,var1,var2,paramsV){
#This function returns the ESS of var, given other variables 
#PIP stands for pairwise invastion plot
#telescoping allows us to zoom in on the ESS by successively decreasing the meshsize and zooming in on the ESS
##############################
		
	#Initialize the meshsize
	minvar = min(V); maxvar = max(V); int = V[2]-V[1]
	INT = (maxvar-minvar)/5
	Vtemp = seq(minvar,maxvar,by=INT)
	
	essV = NULL; teldata = NULL 
	while(INT>int/2){
		outV = try(competeem(var,Vtemp,V,var1,var2,paramsV),silent=TRUE)
		if(is.numeric(outV[[1]])){
			#Save PIP results
			teldata = rbind(teldata,outV$PIPdata) 
			#Reign in infinites
			ess = outV$ess; if(ess==-Inf) ess = minvar; if(ess==Inf) ess = maxvar
			#Zoom in around ess
			INT = INT/2
			Vtemp = seq(ess-INT*2,ess+INT*2,by=INT); 
			if(var=="var1") Vtemp = Vtemp[Vtemp>0]
			if(var=="var2"){Vtemp = Vtemp[Vtemp>=0]; Vtemp = Vtemp[Vtemp<=1]} #this is trait specific!? 
			#Save ESS result
			essV = c(essV,ess)
		}
		#If competeem doesn't work, break the loop. The last working ESS will be the result. 
		if(!is.numeric(outV[[1]])){essV = c(essV,NA); INT = int - .1}
	}
	
	#The final ESS is the last valid result
	essV = essV[!is.na(essV)]; ess = essV[length(essV)]
	
	#Plot the pairwise invasion plot with all of the collected data to monitor working.
	if(plotting) if(!is.null(teldata)) if(is.numeric(teldata[1,1])){
		winners = teldata[teldata[,3]>0,,drop=FALSE]
		#grey is all pairwise sets tested
		plot(teldata[,1],teldata[,2],main=paste(var,"=",ess),xlab="res",ylab="inv",col="gray",pch=20,cex=2)
		#black is where the invasion was successful
		points(winners[,1],winners[,2],pch=20,cex=1.5)
		#green marks the ess found by the algorithm
		points(ess,ess,col="green",pch=20,cex=2)
	}
	return(ess)
}#end PIP_ESS_telescope

##############################
competeem = function(var,Vtemp,V,var1,var2,paramsV){
#function generates a PIP 
##############################

	#setup defaults
	var1.r = var1.i = var1
	var2.r = var2.i = var2
	
	#how far around the resident should the invasions be (I think 1 should always be enough in a deterministic model)
	buff = 1
	
	#matrix for the data; cols: res,inv,invasion potential(>0)
	PIPdata = matrix(ncol=3,nrow=length(Vtemp)*buff*2) 
	inc = Vtemp[2] - Vtemp[1] #the width to space the invaders, used within loop
	i = 0 
	for(res in Vtemp){
		#Set resident value 
		if(var=="var1") var1.r=res; if(var=="var2") var2.r=res
		#Make a Vinv that is centered around res and goes out buff mesh steps but does not include the resident itself. 		
		Vinv = c(seq(res-inc*buff,res-inc,by=inc),seq(res+inc,res+inc*buff,by=inc))
		Vinv = Vinv[Vinv>=min(V)]; Vinv = Vinv[Vinv<=max(V)]
		for(inv in Vinv){
			i = i+1
			#Set the invader value
			if(var=="var1") var1.i=inv; if(var=="var2") var2.i=inv
			Inv = get_inv(var1.r,var1.i,var2.r,var2.i,paramsV=paramsV)
			PIPdata[i,] = c(res,inv,Inv)
		}
	}
	#Take the successful invasions
	winners = PIPdata[PIPdata[,3]>0,,drop=FALSE]; winners = winners[!is.na(winners[,1]),,drop=FALSE]
	#Pick out the ESS from that data
	ess = getESSfromPIPdata(winners)
	
	return(list(PIPdata=PIPdata,ess=ess))  
}#end competeem


##############################
PIP_2D_telescope = function(var1V,var2V,paramsV){
#2D version of PIP_telescope
##############################

	#goal precision levels
	var1int = (var1V[2]-var1V[1])/2; var2int = (var2V[2]-var2V[1])/2
	
	#starting mesh
	starter = 7; var1INT = (max(var1V)-min(var1V))/starter; var2INT = (max(var2V)-min(var2V))/starter
	var1Vstart = seq(min(var1V),max(var1V),by=var1INT); var2Vstart = seq(min(var2V),max(var2V),by=var2INT)
	
	#initalizations
	var1Vgo = var1Vstart; var2Vgo = var2Vstart; minl = 99; ii = 0; reptrack = NULL
	
	outVold = c(999,999)
	outV = c(0,0)
	getout = FALSE
	
	while(minl>1 || var1INT>var1int || var2INT>var2int){
		ii = ii + 1
		outVold = outV
		Mtemp = NULL; var1tick = 0 
		for(var1 in var1Vgo){
			var1tick = var1tick+1; var2tick = 0 
			for(var2 in var2Vgo){
				var2tick = var2tick + 1
				Mtemp = rbind(Mtemp,c(var1,var2,var1tick,var2tick))
			}
		}
		compV = try(competeem2d(Mtemp,min(var1INT,var1int),min(var2INT,var2int),paramsV,plotting=plotting))
		outV.N = compV$outV.N; minl = compV$minl
		if(ii >= 8) minl = 0 #break the loop if neverending
		jitter = 0; if(ii > 4) jitter = runif(1)*.5

		var1INT = var1INT*3/4; var2INT = var2INT*3/4

		if(is.na(outV.N[1,1])){
			var1Vgo = seq(var1Vgo[1]-jitter*var1INT,var1Vgo[length(var1Vgo)]+var1INT,by=var1INT)
			var2Vgo = seq(var2Vgo[1]-jitter*var2INT,var2Vgo[length(var2Vgo)]+var2INT,by=var2INT)
		}
		if(!is.na(outV.N[1,1])){
			outV = c(mean(outV.N[,1]),mean(outV.N[,2]))
			var1Vgo = seq(min(outV.N[,1])-var1INT*5,max(outV.N[,1])+var1INT*5,by=var1INT)
			var2Vgo = seq(min(outV.N[,2])-var2INT*5,max(outV.N[,2])+var2INT*5,by=var2INT)
			outV = c(mean(outV.N[,1]),mean(outV.N[,2]))
			reptrack = rbind(reptrack,outV)
		}
		if(min(var1Vgo)<min(var1V)){var1Vgo=var1Vgo[var1Vgo>=min(var1V)]; var1Vgo = c(min(var1V),var1Vgo)}
		if(max(var1Vgo)>max(var1V)){var1Vgo=var1Vgo[var1Vgo<=max(var1V)]; var1Vgo=c(var1Vgo,max(var1V))}
		if(min(var2Vgo)<min(var2V)){var2Vgo=var2Vgo[var2Vgo>=min(var2V)];var2Vgo=c(min(var2V),var2Vgo)}
		if(max(var2Vgo)>max(var2V)){var2Vgo=var2Vgo[var2Vgo<=max(var2V)];var2Vgo=c(var2Vgo,max(var2V))}
	}
	
	return(list(outV[1],outV[2]))	
}#end PIP_2D_telescope

##############################
competeem2d = function(Mtemp,var1int,var2int,paramsV,plotting=FALSE,presplot=FALSE,xlb="Resident var1",ylb="Resident var2"){
#2D version of competeem
##############################
	
	vectorfield = NULL
	for(res in seq(1,dim(Mtemp)[1])){
		rline = Mtemp[res,,drop=FALSE]
		var1.r = rline[1,1]; var2.r = rline[1,2];
		Minv = matrix(ncol=5,nrow=4)
		Minv[1,] = c(var1.r-var1int,var2.r,-1,0,0)
		Minv[2,] = c(var1.r+var1int,var2.r,1,0,0)
		Minv[3,] = c(var1.r,var2.r-var2int,0,-1,0)
		Minv[4,] = c(var1.r,var2.r+var2int,0,1,0)
		Minv = Minv[Minv[,1]>0,,drop=FALSE]
		Minv = Minv[Minv[,2]>=0,,drop=FALSE]

		for(inv in seq(1,dim(Minv)[1])){
			iline = Minv[inv,,drop=FALSE]
			var1.i = iline[1,1]; var2.i = iline[1,2]
			Inv = get_inv(var1.r,var1.i,var2.r,var2.i,paramsV=paramsV)
			if(!is.na(Inv)) if(Inv>0) Minv[inv,5] = 1
		}
		winners = Minv[Minv[,5]>0,,drop=FALSE]
		if(dim(winners)[1]==0) inv.vector = c(0,0)
		if(dim(winners)[1]>0){
			inv.vector = c(sum(winners[,3]),sum(winners[,4]))
			if(sum(abs(inv.vector)>0)) inv.vector = inv.vector/(inv.vector[1]^2 + inv.vector[2]^2)^.5
		}
		vectorfield = rbind(vectorfield,data.frame(rline,inv.vector[1],inv.vector[2]))		
	}	
	
	if(plotting){
		xint = unique(vectorfield[,1])[2]-unique(vectorfield[,1])[1]
		yint = unique(vectorfield[,2])[2]-unique(vectorfield[,2])[1]		

		xlb = "trait 1"; ylb = "trait 2"
		if(project=="compterr"){xlb=expression(paste("Resident ",italic(r))); ylb=expression(paste("Resident ",italic(x)))}
		plot(vectorfield[,1],vectorfield[,2],type="n",xlab=xlb,ylab=ylb)
		
		vectorfield=vectorfield[!is.na(vectorfield[,5]),]
		vectorField(vectorfield[,5]/yint,vectorfield[,6]/xint,xpos=vectorfield[,1],ypos=vectorfield[,2],scale=0.01,headspan=.1)
		
	}
	
	vectorfield = cbind(vectorfield,seq(1,dim(vectorfield)[1]))
	
	if(dim(vectorfield[vectorfield[,5]^2+vectorfield[,6]^2==0,,drop=FALSE])[1]>0){
		goout = vectorfield[vectorfield[,5]^2+vectorfield[,6]^2==0,,drop=FALSE]
		if(plotting) if(!presplot) points(goout[,1],goout[,2],col="green",pch=20,cex=2)		
		return(list(outV.N=data.frame(col1=goout[,1],col2=goout[,2]),minl=0))
	}
	
	tvar1V = sort(unique(Mtemp[,1]),decreasing=FALSE)
	tvar2V = sort(unique(Mtemp[,2]),decreasing=FALSE)
	tvar1int = tvar1V[2] - tvar1V[1]
	tvar2int = tvar2V[2] - tvar2V[1]

	vdata = NULL; essS = NULL
	for(i in seq(1,dim(vectorfield)[1])){
		line = vectorfield[i,,drop=FALSE]
		line[,1] = line[,1] + tvar1int/2; line[,2] = line[,2] + tvar2int/2
		line[,3] = line[,3] + .5; line[,4] = line[,4] + 0.5
		neighbors = vectorfield[abs(vectorfield[,3]-line[,3])<= .8,]
		neighbors = neighbors[abs(neighbors[,4]-line[,4])<= .8,]
		if(dim(neighbors)[1]==4){
			l = (sum(neighbors[,5])^2 + (sum(neighbors[,6]))^2)^.5
			vdata = rbind(vdata,data.frame(line[,1],line[,2],line[,3],line[,4],l))
		}
		goV = NULL; usum = vsum = 0
	}

	if(min(vdata$l)<3){
		outV.N = vdata[vdata$l==min(vdata$l),,drop=FALSE]
		if(plotting) if(!presplot) points(outV.N[,1],outV.N[,2],col="green",pch=20,cex=2)
		return(list(outV.N=data.frame(col1=outV.N[,1],col2=outV.N[,2]),minl=min(vdata$l)))
	}
	if(min(vdata$l)>3.5){ #everything aligned in one direction
		if(sum(vectorfield[,5])^2 > .5*dim(vectorfield)[1]*2^.5/2)
			if(sum(vectorfield[,6])^2 > .5*dim(vectorfield)[1]*2^.5/2)
				if(sum(vectorfield[,5]<0) && sum(vectorfield[,6]<0)) return(list(outV.N=data.frame(col1=0,col2=0),minl=999))
	}
	if(min(vdata$l)>=2.5) return(list(outV.N=data.frame(col1=NA,col2=NA),minl=999))
}#end competeem2d



##############################
getESSfromPIPdata = function(winners){
#function that returns the ESS given pairwise invasion plot data of only the successful invasions (winners): col1:res, col2: inv
##############################
	amin = NULL #list of the smallest winning invaders that are bigger than the resident
	amax = NULL #list of the biggest winning invaders that are smaller than the resident
	for(res in unique(winners[,1])){
		resdata = winners[winners[,1]==res,,drop=FALSE]
		if(min(resdata[,2])>res) amin = c(amin,min(resdata[,2]))
		if(max(resdata[,2])<res) amax = c(amax,max(resdata[,2]))
	}
	#the ESS is in between the biggest of the smaller winning invaders and the smallest of the bigger winning invaders
	if(is.null(amin)) return(ESS=min(amax,na.rm=TRUE)); if(is.null(amax)) return(ESS=max(amin,na.rm=TRUE))
	
	ESS = mean(c(max(amin,na.rm=TRUE),min(amax,na.rm=TRUE)))
	return(ESS)
}#end getESSfromPIPdata


##############################
test_if_ess = function(var1,var2,var1V,var2V,paramsV,main="",thorough=FALSE){
#here input the var1, var2 values you think are the ESS and see PIPs to check that they are
#gray on the plot indicates that the spot was tested and black indicates that the invader can invade 
#if the ESS you found is correct there should be no black
##############################
	
	var1.r = var1; var2.r = var2
	var1int = var1V[2]-var1V[1] 
	var2int = var2V[2]-var2V[1]
	
	buff = 7	
	
	idata = NULL
	for(var1.i in seq(min(var1V),max(var1V),by=(max(var1V)-min(var1V))/4)){
		for(var2.i in seq(min(var2V),max(var2V),by=(max(var2V)-min(var2V))/4)){
			Inv = get_inv(var1.r,var1.i,var2.r,var2.i,paramsV=paramsV)
			idata = rbind(idata,data.frame(var1.i,var2.i,Inv))
		}
	}
	for(var1.i in seq(var1.r-var1int*buff,var1.r+var1int*buff,by=var1int)){
		for(var2.i in seq(var2.r-var2int*buff,var2.r+var2int*buff,by=var2int)){
			Inv = get_inv(var1.r,var1.i,var2.r,var2.i,paramsV=paramsV)
			idata = rbind(idata,data.frame(var1.i,var2.i,Inv))
		}
	}
	if(thorough){
		for(var1.i in seq(min(var1V),max(var1V),by=(max(var1V)-min(var1V))/100)){
			for(var2.i in seq(min(var2V),max(var2V),by=(max(var2V)-min(var2V))/100)){
				Inv = get_inv(var1.r,var1.i,var2.r,var2.i,paramsV=paramsV)
				idata = rbind(idata,data.frame(var1.i,var2.i,Inv))
			}	
		}
	}
	
	#keeping variables in realistic value, may need changing for different applications
	idata = idata[idata[,1]>0,]
	idata = idata[idata[,2]>=0,]
	
	winner = idata[idata[,3]>0,]
	xlb = "Mutant var1"; ylb = "Mutant var2"
	if(project=="compterr"){ xlb = expression(paste("Mutant ",italic(r))); ylb=expression(paste("Mutant ",italic(x)))}
	plot(idata[,1],idata[,2],xlab=xlb,ylab=ylb,pch=20,col="grey",cex=3,main=main)
	points(winner[,1],winner[,2],pch=20,cex=1.5)
	points(var1.r,var2.r,col="green",pch=20,cex=2)

}#end test_if_ess


##############################
ESSCSSverificationplots_generic = function(essV,var1V,var2V,paramsV){
#function that shows 2D PIPs and vector fields for to allow you to visually verify a solution is both an ESS and a CSS
##############################
		
		X11(width=8,height=3)		
		par(mfrow=c(1,2))

		#plots to show ESS (left panels)
		test_if_ess(essV[[1]],essV[[2]],var1V,var2V,paramsV,main="")

		#plot to show convergence stability
		var1Vgo = seq(min(var1V),max(var1V),length.out=10); var2Vgo = seq(min(var2V),max(var2V),length.out=10)
		Mtemp = NULL; var1tick = 0; 
		for(var1 in var1Vgo){
			var1tick = var1tick+1; var2tick = 0 
			for(var2 in var2Vgo){
				var2tick = var2tick + 1
				Mtemp  = rbind(Mtemp ,c(var1,var2,var1tick,var2tick))
			}
		}
		competeem2d(Mtemp,var1V[2]-var1V[1],var2V[2]-var2V[1],paramsV,plotting=TRUE,presplot=TRUE)
		points(essV[[1]],essV[[2]],col="green",pch=19)

}#end ESSCSSverificationplots_generic
