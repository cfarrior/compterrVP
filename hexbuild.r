#C. Farrior - last updated May 29, 2018 
#cfarrior@gmail.com

#Generic code to build a hexagonal grid and a matrix with information on how close all individuals are to one another if they are on a torus. 

#############################
get_hex = function(size){
#generate two matrices
#pos: a matrix in which each row is an individual (1) x location (2) y location (3) unique number
#neighborsM: a matrix with rows and columns equal to the number of individuals. 
#		each entry gives the shortest distance between the center points of individuals whose unique number is the row number and the individual whose unique number is the column number
#		the grid is wrapped on a torus for these calculations
##############################

	#pos is a matrix of (1) x spatial location (2) y spatial location (3) unique number 
	x.V = seq(1,3*size)
	y.V = seq(1,4*size)

	pos = NULL
	for(y in y.V){
		if(y/2!=floor(y/2)) for(x in seq(2,max(x.V),by=2)) pos = rbind(pos,c(x,y,3)) #odd
		if(y/2==floor(y/2)) for(x in seq(1,max(x.V)-1,by=2)) pos = rbind(pos,c(x,y,4)) #even
	}
	pos[,1] = pos[,1]*cos(pi/6); pos[,2] = pos[,2]*.5 #rescale to make equidistant
	pos[,3] = seq(1,dim(pos)[1]) #unique number for each position
	
	#put this pos on a torus and give the neighbor distances
	#neighborsM is a big, symmetric matrix with all individuals as the rows and columns. 
	#entries will be the distance between the two individuals when that individual is moved to the center of the hexagonal grid. 	
	maxx = max(pos[,1])+0.1; minx = min(pos[,1])
	maxy = max(pos[,2])+0.1; miny = min(pos[,2])
	center = closest(mean(pos[,1]),mean(pos[,2]),pos)
	
	neighborsM = matrix(nrow=dim(pos)[1],ncol=dim(pos)[1])
	for(i in seq(1,dim(pos)[1])){
		line = pos[i,,drop=FALSE]
		
		temp = pos
		temp[,1] = temp[,1] - line[,1] + center[,1]
		temp[,2] = temp[,2] - line[,2] + center[,2]

		#now we need to torus 
		#(1) move individuals that sitck out on the right to the left
		overx = temp[temp[,1]>maxx,]; if(dim(overx)[1]>0){
			temp = temp[temp[,1]<=maxx,]
			overx[,1] = overx[,1]-min(overx[,1])+minx
			temp=rbind(temp,overx)
		}
		#(2) move individuals that sitck out on the top to the bottom
		overy = temp[temp[,2]>maxy,]; if(dim(overy)[1]>0){
			temp = temp[temp[,2]<=maxy,]
			overy[,2] = overy[,2]-min(overy[,2])+miny
			temp = rbind(temp,overy)
		}
		#(3) move individuals that sitck out on the left to the right
		underx = temp[temp[,1]<minx,]; if(dim(underx)[1]>0){
			temp = temp[temp[,1]>=minx,]
			underx[,1] = underx[,1]-min(underx[,1])+max(temp[,1])+cos(pi/6)
			temp = rbind(temp,underx)
		}
		#(4) move individuals that sitck out on the bottom to the top
		undery = temp[temp[,2]<miny,]; if(dim(undery)[1]>0){
			temp = temp[temp[,2]>=miny,]
			undery[,2] = undery[,2]-min(undery[,2])+max(temp[,2])+.5
			temp = rbind(temp,undery)
		}
			
		for(j in seq(1,dim(temp)[1])){
			jline = temp[j,,drop=FALSE]
			neighborsM[line[,3],jline[,3]] = ((jline[,1]-center[,1])^2+(jline[,2]-center[,2])^2)^.5
		}
	}	
	
	return(list(pos=pos,neighborsM=neighborsM))	
}


##############################
closest = function(x,y,tpos=NaN){
#returns the individual closest to the x, y coordinates given
#this only works if you don't need to torus
##############################
	if(is.na(tpos[1])) tpos = pos
	distV = (tpos[,1]-x)^2 + (tpos[,2]-y)^2
	return(tpos[distV==min(distV),,drop=FALSE][1,,drop=FALSE])
}#end closest


##############################
get_neighbors = function(line,data,dist){
#return all neighbors of a target individual closer than a distance, dist
#		excludes the individual itself
#line: the target individual
#data: the positions of all individuals
#dist: the distance threshold
##############################
	
	neighdata = cbind(data,neighborsM[line[,3],])
	
	neighdata = neighdata[neighdata[,dim(neighdata)[2]]<dist,,drop=FALSE]
	neighdata = neighdata[neighdata[,3]!=line[,3],,drop=FALSE]
	
	return(neighdata)	
}#end get_neighbors


##############################
get_neighbors_traditional = function(line,data,dist){
##############################
	
	x = line[1]; y = line[2] 
	
	neighdata = cbind(data,NaN); distcol = dim(neighdata)[2]
	neighdata[,distcol] = (neighdata[,1] - x)^2 + (neighdata[,2]-y)^2
	neighdata = neighdata[neighdata[,distcol]<dist^2,,drop=FALSE]
	neighdata = neighdata[neighdata[,3]!=line[3],,drop=FALSE]
	
	return(neighdata)
	
}#end get_neighbors_traditional
