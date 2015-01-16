#' Enumeration of all possible tables
#'
#' A Function to enumerate all tables for N decisions in a special form
#' @param N is integer for number of patients to enter.
#' @param X is four values of 2x2 table to start from. Default is c(0,0,0,0)
#' @return Xsm is a list of length N+1, each of which is a list of matrices. The values of matrices are 0 throughout.
#' @return Ysm is a list of matrix, each row of which is 4 values of 2x2 table.
#' @return M is an utility output.
#' @return Z is a list of 4 elements for Aarm&Success, Aarm&Failure, Barm&Succsess, and Barm&Failure; each of which is in the same form with X.
#' @keywords NA
#' @export
#' @examples
#' N <- 20
#' out <- serialTables(N)


#library(combinat)
# args N : number of max patients
# values Xsm : Container of value 0 of all tables of 0,1,2,...,N
# values Ysm : 4 cell values of all tables of 0,1,2,...,N
serialTables <- function(N,X0=rep(0,4)){
  # Number of patients allocated to each arm
  M <- list()
  for(i in 0:N){
    M[[i+1]] <- as.matrix(xsimplex(2,i)) # xsimplex() is to generate all points on a(p,n) simplex lattice.
  }
  # Total sum of probability of three types of decision; armA, armB and by-chance-AorB
  decision <- rep(0,3)
  
  # Xsm is the list of probability of all tables
  # Xsm[i] is when number of total patients is i-1
  # It is consisted of (1*(i+1)),(2*i),(3*(i-1)),...,(i*2),((i+1)*1) matrices
  # For each table, 4 elements should be assigned as an array Ysm.
  Xsm <- Ysm <- ZAS <- ZAF <- ZBS <- ZBF <- list()
  ZAS[[1]] <- matrix(X0[1],1,1)
  ZAF[[1]] <- matrix(X0[2],1,1)
  ZBS[[1]] <- matrix(X0[3],1,1)
  ZBF[[1]] <- matrix(X0[4],1,1)
  Xsm[[1]] <- list(matrix(0,1,1))
  Ysm[[1]] <- list(matrix(rep(0,4),nrow=1))
  for(i in 1:N){
    n <- i+1
    Xsm[[n]] <- list()
    ZAS[[n]] <- ZAF[[n]] <- ZBS[[n]] <- ZBF[[n]] <- list()
    for(j in 1:n){
      Xsm[[n]][[j]] <- matrix(0,j,n-j+1)
      ZAS[[n]][[j]] <- matrix(rep((j-1):0,n-j+1),j,n-j+1)
      ZAF[[n]][[j]] <- (j-1) - ZAS[[n]][[j]]
    }
    for(j in 1:n){
      ZBS[[n]][[j]] <- t(ZAS[[n]][[n+1-j]])
      ZBF[[n]][[j]] <- t(ZAF[[n]][[n+1-j]])
		}
    Ysm[[n]] <- list()
    for(j in 1:n){
      Y1 <- t(M[[j]])
      Y2 <- t(M[[n-j+1]])
      cmb <- expand.grid(1:(j),1:(n-j+1))
      tmp <- cbind(matrix(Y1[cmb[,1],],ncol=2),matrix(Y2[cmb[,2],],ncol=2))
      Ysm[[n]][[j]] <- t(t(tmp)+X0)
    }
  }
  Z <- list(ZAS,ZAF,ZBS,ZBF) # A&Success,A&Failure,B&Success, B&Failure counts in the same data structure with Xsm
  return(list(Xsm,Ysm,M,Z))
}

#' Utility function
#'
#' Utility
#' @param m a vector of 4 elements
#' @return success rate
#' @keywords NA
#' @export
#' @examples
#' calc.success2(c(1,2,3,4))

calc.success2 <- function(m){
	fl <- m[,2]+m[,4]
	suc <- m[,1]+m[,3]
	suc/(suc+fl)
}

#' Utility function2
#'
#' Utility2
#' @param m a vector of 4 elements
#' @return A-selection rate
#' @keywords NA
#' @export
#' @examples
#' calc.select2(c(1,2,3,4))

calc.select2 <- function(m){ # sum(m)!=0
	a <- m[,1]+m[,2]
	b <- m[,3]+m[,4]
	a/(a+b)
}


#' Expected value of beta-distribution-dependent decision 
#'
#' Expected value of beta-distribution-dependent decision 
#' @param x a vector of 4 elements
#' @return probability to select A 
#' @keywords NA
#' @export
#' @examples
#' prob.expected.3(c(1,2,3,4))
# Decision rule gives prob to select arm A based on 4 numbers of each table.

# Calculate probability for conventional (expected value-based) to select A 

prob.expected.3 <- function(x){ 
	x. <- x+1
	a <- x.[,1]/(x.[,1]+x.[,2])
	b <- x.[,3]/(x.[,3]+x.[,4])
	(sign(a-b)+1)/2
   
}


#' Probability of better than target success rate-dependent self decision 
#'
#' Probability of better than target success rate-dependent self decision 
#' @param x a vector of 4 elements
#' @return probability to select A 
#' @keywords NA
#' @export
#' @examples
#' prob.self.3(c(1,2,3,4))

# Decision rule gives prob to select arm A based on 4 numbers of each table.

# Calculate probability for targetting people to select A 
prob.self.3 <- function(x,target){ 
	x. <- x+1
	a <- pbeta(target,x.[,1],x.[,2],lower.tail=FALSE)
	b <- pbeta(target,x.[,3],x.[,4],lower.tail=FALSE)
	(sign(a-b)+1)/2
}

#' Make a list of selection probability
#'
#' Make a list of selection probability
#' @param tablist list of matrices with 4 columns
#' @param f function to calculate decision probability
#' @return list of vectors
#' @keywords NA
#' @export
#' @examples
#' N <- 20
#' out <- serialTables(N)
#' prob.success.v <- prob.list(out[[2]],calc.success2)
#' prob.select.A <- prob.list(out[[2]],calc.select2)

prob.list <- function(tablist,f){
	prob.success.v <- list()
	for(i in 1:length(tablist)){
		prob.success.v[[i]] <- lapply(tablist[[i]],f)
	}
	prob.success.v[[1]][[1]] <- 0.5
	prob.success.v
}

#' Utility function to handle lists
#'
#' Utility function to handle lists
#' @param L1 list
#' @param L2 list with the same form as L1
#' @return list 
#' @keywords NA
#' @export
#' @examples
#' L1 <- list(c(1,2,3),matrix(1:4,2,2))
#' L2 <- list(c(2,3,4),matrix(2:5,2,2))
#' relist.sum(L1,L2)

relist.sum <- function(L1,L2){
  L1.a <- as.relistable(L1)
  L2.a <- as.relistable(L2)
  L1.v <- unlist(L1.a)
  L2.v <- unlist(L2.a)
  ret <- L1.v + L2.v
  return(relist(ret))
}

#' From the output of serialTables(), multiple lists are generated
#'
#' From the output of serialTables(), multiple lists are generated
#' @param out; serialTables() ' out
#' @return list of success rate for tables
#' @return list of selection A rate for tables
#' @return list of A&Success, list of A&Failure, list of B&Success and list of B&Failure
#' @keywords NA
#' @export
#' @examples
#' N <- 20
#' out <- serialTables(N)
#' outseries <- unlistSerialTables(out)

unlistSerialTables <- function(out){
	AS <- lapply(out[[4]][[1]],unlist)
	AF <- lapply(out[[4]][[2]],unlist)
	BS <- lapply(out[[4]][[3]],unlist)
	BF<- lapply(out[[4]][[4]],unlist)
	prob.success.v <- prob.list(out[[2]],calc.success2)
	prob.select.A <- prob.list(out[[2]],calc.select2)
	success.series <- lapply(prob.success.v,unlist)
	selectA.series <- lapply(prob.select.A,unlist)
	# 2D grid of success and selection
	success.select.series <- list()
	N <- length(success.series)
	for(i in 1:N){
		success.select.series[[i]] <- round(cbind(success.series[[i]],selectA.series[[i]])*(i-1))
	}

	return(list(AS=AS,AF=AF,BS=BS,BF=BF,success.series=success.series,selectA.series=selectA.series,ssgrid = success.select.series))
}


#' Make a list of expect-dependent selection probability
#'
#' Make a list of expect-dependent selection probability
#' @param tablist list of matrices with 4 columns
#' @return list of selection probability
#' @keywords NA
#' @export
#' @examples
#' N <- 20
#' out <- serialTables(N)
#' prob.select.exp <- prob.expect.select(out[[2]])

prob.expect.select <- function(tablist){
	n <- length(tablist)
	prob.select.exp <- list()
	for(i in 1:n){
		tmp <- lapply(tablist[[i]],prob.expected.3)
		prob.select.exp[[i]] <- list()
		for(j in 1:length(tmp)){
			prob.select.exp[[i]][[j]] <- matrix(tmp[[j]],j,i-j+1)
		}
	}
	prob.select.exp
}

#' Make a list of self-decision selection probability
#'
#' Make a list of self-decision selection probability
#' @param tablist list of matrices with 4 columns
#' @param target a value from 0 to 1 
#' @return list of selection probability
#' @keywords NA
#' @export
#' @examples
#' N <- 20
#' out <- serialTables(N)
#' targets <- seq(from=0,to=1,by=0.1)
#' prob.select.selfs <- list()
#' for(i in 1:length(targets)){
#'  prob.select.selfs[[i]] <- prob.self.select(out[[2]],targets[i])
#' }

prob.self.select <- function(tablist,target){
	n <- length(tablist)
	prob.select.exp <- list()
	for(i in 1:n){
		tmp <- lapply(tablist[[i]],prob.self.3,target)
		prob.select.exp[[i]] <- list()
		for(j in 1:length(tmp)){
			prob.select.exp[[i]][[j]] <- matrix(tmp[[j]],j,i-j+1)
		}
	}
	prob.select.exp
}

#' Make a list of selection probability for a mixture of selection types
#'
#' Make a list of selection probability for a mixture of selection types
#' @param list of selection probabilities
#' @param a weight vector
#' @return list of selection probability
#' @keywords NA
#' @export
#' @examples
#' N <- 20
#' out <- serialTables(N)
#' prob.select.exp <- prob.expect.select(out[[2]])
#' targets <- seq(from=0,to=1,by=0.1)
#' prob.select.selfs <- list()
#' for(i in 1:length(targets)){
#'  prob.select.selfs[[i]] <- prob.self.select(out[[2]],targets[i])
#' }
#' prob.selection.list <- prob.select.selfs
#' prob.selection.list[[length(prob.selection.list)+1]] <- prob.select.exp
#' weight.v <- runif(length(prob.selection.list))
#' weight.v <- weight.v/sum(weight.v)
#' prob.selection.this <- prob.select.mix(prob.selection.list,weight.v)

prob.select.mix <- function(L,v){
	prob.selection.this <- list()
	N <- length(L[[1]])
	for(i in 1:N){
		prob.selection.this[[i]] <- relist(unlist(as.relistable(L[[1]][[i]]))*v[1])
		for(j in 2:length(v)){
			tmp <- relist(unlist(as.relistable(L[[j]][[i]]))*v[j])
			prob.selection.this[[i]] <- relist.sum(prob.selection.this[[i]], tmp)
		}
	}
	prob.selection.this
}

#' Probability of all tables for particular selection pattern and success rate fo A and B arms
#'
#' Probability of all tables for particular selection pattern and success rate fo A and B arms
#' @param list of A-selection probability
#' @param a vector of two arms success rate
#' @return list of list of probability of tables
#' @keywords NA
#' @export
#' @examples
#' #### Step 1
#' N <- 20
#' out <- serialTables(N)
#' outseries <- unlistSerialTables(out)
#'
#' #### Step 2
#' prob.select.exp <- prob.expect.select(out[[2]])
#' targets <- seq(from=0,to=1,by=0.1)
#' prob.select.selfs <- list()
#' for(i in 1:length(targets)){
#'  prob.select.selfs[[i]] <- prob.self.select(out[[2]],targets[i])
#' }
#' prob.select.list <- prob.select.selfs
#' prob.select.list[[length(prob.select.list)+1]] <- prob.select.exp
#' 
#' #### Step 3
#' pop1 <- c(0,0,0,0,0,0,0,0,0,0,0,1)
#' select1 <- prob.select.mix(prob.select.list,pop1)
#' p.AB <- c(0.8,0.5)
#' tab.prob <- table.prob(select1,p.AB)

table.prob <- function(s1,p){
	p.A <- p[1]
	p.B <- p[2]
	ret <- s1
	ret[[1]][[1]] <- 1
	N <- length(s1)-1
	for(i in 1:N){
		n <- i+1
		ret[[n]] <- lapply(ret[[n]],"*",0)
		for(j in 1:i){
			tmp.select <- s1[[i]][[j]]
			dims <- dim(tmp.select)
			# When A is selected j+1, because the submatrix-number is based on the number of A-selections.
			# When A.Failure, row number increases
			# When B.Failure, column number increases
			# (1) A & Success
			J <- j+1;xspan <- 1:dims[1];yspan <- 1:dims[2];
			ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + tmp.select * ret[[i]][[j]] * p.A
			# (2) A & Failure
			J <- j+1;xspan <- 2:(dims[1]+1);yspan <- 1:dims[2];
			ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + tmp.select * ret[[i]][[j]] * (1-p.A)
			# (3) B & Success
			J <- j;xspan <- 1:dims[1];yspan <- 1:dims[2];
			ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + (1-tmp.select) * ret[[i]][[j]]* p.B
			# (4) B & Failure
			J <- j;xspan <- 1:dims[1];yspan <- 2:(dims[2]+1);
			ret[[n]][[J]][xspan,yspan] <- ret[[n]][[J]][xspan,yspan] + (1-tmp.select) * ret[[i]][[j]] * (1-p.B)
		}
	}
	ret
}


#' Weighted mean for lists
#'
#' Weighted mean for lists
#' @param v list of value vector
#' @param w list of weight vector
#' @return vector of weighted mean
#' @keywords NA
#' @export
#' @examples
#' #### Step 1
#' N <- 20
#' out <- serialTables(N)
#' outseries <- unlistSerialTables(out)
#'
#' #### Step 2
#' prob.select.exp <- prob.expect.select(out[[2]])
#' targets <- seq(from=0,to=1,by=0.1)
#' prob.select.selfs <- list()
#' for(i in 1:length(targets)){
#'  prob.select.selfs[[i]] <- prob.self.select(out[[2]],targets[i])
#' }
#' prob.select.list <- prob.select.selfs
#' prob.select.list[[length(prob.select.list)+1]] <- prob.select.exp
#' 
#' #### Step 3
#' pop1 <- c(0,0,0,0,0,0,0,0,0,0,0,1)
#' select1 <- prob.select.mix(prob.select.list,pop1)
#' p.AB <- c(0.8,0.5)
#' tab.prob <- table.prob(select1,p.AB)
#' tab.prob.unlisted <- lapply(tab.prob,unlist)
#' mean.successes <- weighted.mean.list(outseries$success.series,tab.prob.unlisted)

weighted.mean.list <- function(v,w){
	n <- length(v)
	ret <- rep(0,n)
	for(i in 1:n){
		ret[i] <- weighted.mean(v[[i]],w[[i]])
	}
	ret
}

#' mean of success rate and selection A
#'
#' mean of success rate and selection A
#' @param outerseries output
#' @param tab.prob.unlisted
#' @return a matrix consisted of 2 column vectors
#' @keywords NA
#' @export
#' @examples
#' #### Step 1
#' N <- 20
#' out <- serialTables(N)
#' outseries <- unlistSerialTables(out)
#'
#' #### Step 2
#' prob.select.exp <- prob.expect.select(out[[2]])
#' targets <- seq(from=0,to=1,by=0.1)
#' prob.select.selfs <- list()
#' for(i in 1:length(targets)){
#'  prob.select.selfs[[i]] <- prob.self.select(out[[2]],targets[i])
#' }
#' prob.select.list <- prob.select.selfs
#' prob.select.list[[length(prob.select.list)+1]] <- prob.select.exp
#' 
#' #### Step 3
#' pop1 <- c(0,0,0,0,0,0,0,0,0,0,0,1)
#' select1 <- prob.select.mix(prob.select.list,pop1)
#' p.AB <- c(0.8,0.5)
#' tab.prob <- table.prob(select1,p.AB)
#' tab.prob.unlisted <- lapply(tab.prob,unlist)
#' mean.success.select <- mean.success.select(outseries,tab.prob.unlisted)
#' matplot(mean.success.select,type="l")

mean.success.select <- function(os,w){
	n <- length(w)
	ret1 <- ret2 <- rep(0,n)
	for(i in 1:n){
		ret1[i] <- weighted.mean(os$success.series[[i]],w[[i]])
		ret2[i] <- weighted.mean(os$selectA.series[[i]],w[[i]])
	}
	return(cbind(ret1,ret2))
}


#' 2D-grid probability of Success rate, SelectionA rate
#'
#' 2D-grid probability of Success rate, SelectionA rate
#' @param outerseries output
#' @param tab.prob.unlisted
#' @return a list of matrix (x axis is success rate and y axis is selection.A rate)
#' @keywords NA
#' @export
#' @examples
#' #### Step 1
#' N <- 20
#' out <- serialTables(N)
#' outseries <- unlistSerialTables(out)
#'
#' #### Step 2
#' prob.select.exp <- prob.expect.select(out[[2]])
#' targets <- seq(from=0,to=1,by=0.1)
#' prob.select.selfs <- list()
#' for(i in 1:length(targets)){
#'  prob.select.selfs[[i]] <- prob.self.select(out[[2]],targets[i])
#' }
#' prob.select.list <- prob.select.selfs
#' prob.select.list[[length(prob.select.list)+1]] <- prob.select.exp
#' 
#' #### Step 3
#' pop1 <- c(0,0,0,0,0,0,0,0,0,0,0,1)
#' select1 <- prob.select.mix(prob.select.list,pop1)
#' p.AB <- c(0.8,0.5)
#' tab.prob <- table.prob(select1,p.AB)
#' tab.prob.unlisted <- lapply(tab.prob,unlist)
#' dist.ss.out <- dist.success.select(outseries,tab.prob.unlisted)
#' for(i in 1:length(dist.ss.out)){
#' 	persp(dist.ss.out[[2]][[i]],theta=80,phi=50,shade=0.75)
#' }

dist.success.select <- function(os,prob.v){
	ss.grid <- os$ssgrid
	N <- length(prob.v)
	grid.prob.series <- list()
	for(i in 1:N){
		grid.prob.series[[i]] <- matrix(0,i,i)
		for(j in 1:length(ss.grid[[i]][,1])){
			grid.prob.series[[i]][ss.grid[[i]][j,1]+1,ss.grid[[i]][j,2]+1] <- grid.prob.series[[i]][ss.grid[[i]][j,1]+1,ss.grid[[i]][j,2]+1] + prob.v[[i]][j]
		}
	}
	return(list(ss.grid=ss.grid,pr.grid=grid.prob.series))
}

