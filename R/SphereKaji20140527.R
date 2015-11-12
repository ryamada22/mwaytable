##########
### Handlings of multi-way tables
##########
# Marginal count
# A is an array
#' Calculate marginal counts of array
#'
#' @param A An array or matrix
#' @return A list of marginal counts per variable
#' @examples
#' r <- c(2,3,4)
#' A <- array(1:prod(r),r)
#' calc.marg(A)
calc.marg <- function(A){
  d <- dim(A)
  ret <- list()
  for(i in 1:length(d)){
    ret[[i]] <- apply(A,i,sum)
  }
  ret
}
# Expected table
#' Make an expected table
#'
#' @param A an array or matrix
#' @return The expected table in a shape of array or matrix
#' @examples
#' r <- c(2,3,4)
#' A <- array(1:prod(r),r)
#' make.exp.table(A)
make.exp.table <- function(A){
  n <- sum(A)
  marg <- calc.marg(A)
  tmp <- marg[[1]]
  for(i in 2:length(marg)){
    tmp <- t(marg[[i]]/n) %x% tmp
  }
  tmp
}
# Difference table
#' Make a table, the given table minus its expected table
#'
#' @param A an array or matrix
#' @return The difference table in a shape of array or matrix
#' @examples
#' r <- c(2,3,4)
#' A <- array(1:prod(r),r)
#' make.diff.table(A)
make.diff.table <- function(A){
  E <- make.exp.table(A)
  array(c(A)-c(E),dim(A))
}
# Chisquare value
#' Calculate chi-square value of an array or matrix
#'
#' @param A an array or matrix
#' @return The chi-square value of the table
#' @examples
#' r <- c(2,3,4)
#' A <- array(1:prod(r),r)
#' calc.chisq(A)
calc.chisq <- function(A){
  E <- make.exp.table(A)
  D <- make.diff.table(A)
  sum(c(D)^2/c(E))
}
# Simplex and rotation of simplex
# (d-1)-simplex is an polytope with d vertices in d-dimensional space
# This function gives coordinates of d vertices in d-1 dimensional space
#' Coordinate of simplex vertices
#'
#' @param d Number of vertices of a simplex
#' @return d times (d-1) matrix, each row of which is a coordinate of a vertex of (d-1)-simplex, whose center is the origin and distance of whose vertices from the origin is 1
#' @examples
#' d <- 3
#' cv <- CategoryVector(d)
#' print(cv)
#' apply(cv^2,1,sum)
#' cv %*% t(cv) 
CategoryVector<-function(d){
  df <- d - 1
  diagval <- 1:d
  diagval <- sqrt((d)/df) * sqrt((d - diagval)/(d - diagval + 1))
  others <- -diagval/(d - (1:d))
  m <- matrix(rep(others, d), nrow = d, byrow = TRUE)
  diag(m) <- diagval
  m[upper.tri(m)] <- 0
  as.matrix(m[, 1:df])
}

# This function gives coordinates of d vertices in "d" dimensional space with d-th coordinate of all vertices being equal and distance of all vertices from origin being 1.
# This is a rotation matrix.
# We call the output "simplex-rotation matrix".
#' Make a rotation matrix using simple vertices coordinates (1)
#'
#' @param k A integer
#' @return k x k rotation matrix to transfer k vertices so that all k vertices' k-th coordinates are the same
#' @examples
#' k <- 3
#' M <- make.simplex.0(k)
#' M %*% diag(rep(1,k))
make.simplex.0 <- function(k){
  cv <- CategoryVector(k)
  rbind(t(cv*sqrt(1-1/k)),rep(1/sqrt(k),k))
}

#' Make a rotation matrix using simple vertices coordinates (2)
#'
#' @param k A integer
#' @return k x k rotation matrix to transfer k vertices so that all k vertices' k-th coordinates are the same
#' @examples
#' k <- 3
#' M <- make.simplex(k)
#' M %*% diag(rep(1,k))
make.simplex <- function(k){
  ret <- matrix(0,k,k)
  for(i in 1:(k-1)){
    for(j in 1:k){
      if(j < i){
      }else if(j==i){
        ret[i,j] <-  sqrt((k-i)/(k-i+1))
      }else{
        ret[i,j] <- -sqrt(1/((k-i)*(k-i+1)))
      }
    }
  }
  ret[k,] <- sqrt(1/k)
  ret
}

# Kronecker product of multiple simplex-rotation matrices
#' Kronecker product of simplex-based rotation matrices 
#'
#' @param r A vector of positive integers, a vector of number of levels of variables
#' @return Kronecker product of rotation matrices, each of which is a simplex-based rotation matrix of number of levels of each variable
#' @examples
#' r <- c(2,3,4)
#' KM <- make.simplex.multi(r)
#' dim(KM)
make.simplex.multi <- function(r){
  X <- make.simplex(r[1])
  k <- length(r)
  if(k > 1){
    for(i in 2:k){
      X <- make.simplex(r[i]) %x% X
    }
  }
  X
}
# This function's output indicates which elements of table vectors should be zero after the kronecker product rotation.
# r is dimension vector or multi-way table
#' Make a vector indexing rows that should be zero
#'
#' @param r A vector of positive integers, a vector of number of levels of variables
#' @return A vector with 0 or 1; indices with 1 means that elements of table vectors should be zero after the kronecker product rotation.
#' @examples
#' r <- c(2,3,4)
#' v0 <- make.vector0(r)

make.vector0 <- function(r){
  tmp <- list()
  for(i in 1:length(r)){
    tmp[[i]] <- c(rep(1,r[i]-1),0)
  }
  tmp2 <- apply(expand.grid(tmp),1,sum)
  as.numeric(tmp2>1)
}

# r is dimension vector or multi-way table
# Kronecker-product rotation matrix: X
# The meaningfull rows of X: X.sub
# The X^(-1) = X^T: X.inv
# The X.sub^(-1): X.inv.sub
# The vector indicating which rows are meaningfull/less
#' Multiway table-related rotation matrix and its features
#'
#' @param r A vector of positive integers, a vector of number of levels of variables
#' @return X A Rotation matrix; output of make.simplex.multi()
#' @return X.sub A matrix only with meaningful rows of X
#' @return X.inv Inverse matrix of X
#' @return X.inv.sub Meaningful part of X.inv
#' @return vector0 A vector indicating meaningful rows of X
#' @examples
#' r <- c(2,3,4)
#' mXv <- make.X.vector0(r)
#' mXv$X %*% mXv$X.inv
make.X.vector0 <- function(r){
  X <- make.simplex.multi(r)
  #X.inv <- solve(X)
  X.inv <- t(X)
  vector0 <- make.vector0(r)
  non.zero <- which(vector0==1)
  return(list(X=X,X.sub=X[non.zero,],X.inv=X.inv,X.inv.sub=X.inv[,non.zero],vector0=vector0))
}
###
# multi-way table is vectorized into a R-length vector and it is rotated with map2df matrix to a vector whose R-df elements are zero.
#' Make a matrix that transfer a vectorized table array to a vector in the df-dimensional space
#'
#' @param A An array or matrix of multiway table
#' @return A matrix
#' @examples
#' r <- c(2,3,4)
#' A <- array(1:prod(r),r)
#' M <- map2df(A)
#' M %*% c(A)
map2df <- function(A){
  r <- dim(A)
  out <- make.X.vector0(r)
  E.A <- make.exp.table(A)
  W <- t(out$X.inv.sub) %*% diag(1/c(E.A)) %*% (out$X.inv.sub)
  eigen.W <- eigen(W)
  diag(sqrt(eigen.W[[1]])) %*% t(eigen.W[[2]]) %*% out$X.sub
}
#' Make a matrix that transfer back a vector in the df-dimensional space to a vector with number-of-cell of table elements
#'
#' @param A An array or matrix of multiway table
#' @return A matrix
#' @examples
#' r <- c(2,3,4)
#' A <- array(1:prod(r),r)
#' M <- map2df(A)
#' M.inv <- map2full(A)
#' M.inv %*% (M %*% c(A))
map2full <- function(A){
  r <- dim(A)
  out <- make.X.vector0(r)
  E.A <- make.exp.table(A)
  W <- t(out$X.inv.sub) %*% diag(1/c(E.A)) %*% (out$X.inv.sub)
  eigen.W <- eigen(W)
  out$X.inv.sub %*% eigen.W[[2]] %*% diag(1/sqrt(eigen.W[[1]]))
}
# A test vector in R-dimensional space is roted to df-dimensional vector with this matrix
#' Make a matrix that transfer test vectors to the direction in the df-dimensional space
#'
#' @param A An array or matrix of multiway table
#' @return A matrix
#' @examples
#' r <- c(2,3,4)
#' A <- array(1:prod(r),r)
#' M <- map2dfdir(A)
#' c(A) %*% M
map2dfdir <- function(A){
  r <- dim(A)
  out <- make.X.vector0(r)
  E.A <- make.exp.table(A)
  W <- t(out$X.inv.sub) %*% diag(1/c(E.A)) %*% (out$X.inv.sub)
  eigen.W <- eigen(W)
  out$X.inv.sub %*% (eigen.W[[2]]) %*% diag(1/sqrt(eigen.W[[1]]))
}
###
# 3 matrices mentioned above are bundled.
#' Make three matrices that are the output of map2df(), map2full() and map2dfdir
#'
#' @param A An array or matrix of multiway table
#' @return Full2DF A matrix transferring table vectors to df-dimension vectors
#' @return DF2full A matrix transferring back vectors in df-dimensional space back to vectors with number-of-cell elements
#' @return Normalvec2DF A matric transferring test vectors to df-dimension vectors
#' @examples
#' r <- c(2,3,4)
#' A <- array(1:prod(r),r)
#' map.matrix(A)

map.matrix <- function(A){
  r <- dim(A)
  out <- make.X.vector0(r)
  E.A <- make.exp.table(A)
  W <- t(out$X.inv.sub) %*% diag(1/c(E.A)) %*% (out$X.inv.sub)
  eigen.W <- eigen(W)
  Full2DF <- diag(sqrt(eigen.W[[1]])) %*% t(eigen.W[[2]]) %*% out$X.sub
  DF2full <- out$X.inv.sub %*% eigen.W[[2]] %*% diag(1/sqrt(eigen.W[[1]]))
  NormalVec2DF <- out$X.inv.sub %*% (eigen.W[[2]]) %*% diag(1/sqrt(eigen.W[[1]]))
  return(list(Full2DF = Full2DF,DF2full = DF2full, NormalVec2DF = NormalVec2DF))
}
# multi-dimensional random vectors in normal distribution whose length is 1. 
#' Generate random unit vectors
#'
#' @param n Number of random vecotrs to be generated
#' @param df An integer indicating dimension
#' @return A n x df matrix, each row of which is a random unit vector
#' @examples
#' n <- 100
#' df <- 2
#' plot(st.mvn(n,df))
st.mvn <- function(n,df){
  R <- matrix(rnorm(n*df),ncol=df)
  L <- sqrt(apply(R^2,1,sum))
  R/L
}


# integrate segments with "OR"
#' Union of segment(s) of multiple one-dimensional segments
#'
#' @param s A matrix with two columns each row of which indicating start and end of each segment
#' @return A matrix with two columns each row of which indicating unified segment
#' @examples
#' s <- matrix(c(0,4,1,6,2,3),byrow=TRUE,ncol=2)
#' segment.union(s)

segment.union <- function(s){
  n.seg <- length(s[,1])
  my.seg.1 <- t(apply(s,1,sort))
  ord <- order(my.seg.1[,1])
  my.seg.2 <- my.seg.1[ord,]
  st.end <- rep(0:1,n.seg)
  st.end
  loc <- c(t(my.seg.2))
  loc.ord <- order(loc)
  loc.ordered <- loc[loc.ord]
  st.end.ordered <- st.end[loc.ord]
  
  new.st <- c()
  current <- 0
  cnt <- 0
  ret <- matrix(0,0,2)
  for(i in 1:length(loc)){
    if(st.end.ordered[i]==0){
      if(length(new.st)==0){
        new.st <- loc.ordered[i]
      }
      current <- current + 1
    }else{
      if(current==1){
        cnt <- cnt+1
        #ret[[cnt]] <- c(new.st,loc.ordered[i])
        ret <- rbind(ret,c(new.st,loc.ordered[i]))
        new.st <- c()
      }
      current <- current - 1
    }
  }
  ret
}
# integrate segments with "AND"
#' Overlapped segment(s) of multiple one-dimensional segments
#'
#' @param s A matrix with two columns each row of which indicating start and end of each segment
#' @return A matrix with two columns each row of which indicating overlapped segment
#' @examples
#' s <- matrix(c(0,4,1,6,2,3),byrow=TRUE,ncol=2)
#' segment.overlap(s)
segment.overlap <- function(s){
  n.seg <- length(s[,1])
  my.seg <- t(apply(s,1,sort))
  max1 <- max(my.seg[,1])
  min2 <- min(my.seg[,2])
  ret <- c(max1,min2)
  if(max1>min2){
    ret <- c(0,0)
  }
  ret
}

# Chisquare distribution for +chi^2 and (-chi)^2 and cumulative function from -inf -> 0 -> +inf
#' Cumulative distribtution of bidirectional chisquare distribution supporting -Inf to +Inf 
#'
#' @param a Chi-square value or its negative
#' @param df Degrees of freedom
#' @return Cumulative probability
#' @examples
#' as <- seq(from=-5,to=5,length=101)
#' plot(pchisq.bid(as,1))
pchisq.bid <- function(a,df){
  tmp <- as.numeric(a>0) - sign(a)*pchisq(a^2,df,lower.tail = FALSE)/2
  tmp[which(a==0)] <- 0.5
  tmp
}

###
# Marginal counts-dependent matarials
# Arguments
#  A multi-way table (array)
# Values
#  tables=list(data=A,expected=exp.table,diff=diff.table)
##  A: input array, expected: expected array, diff: diff array
#  r.vec = r.vec
## Dimensions of array
#  df=df
## degree of freedom
#  matrices = list(X=out$X,X.sub=out$X.sub,X.inv=out$X.inv,X.inv.sub=out$X.inv.sub)
## matrices to rotate simplex between full and df dimensional spaces
#  zero.vec=out$vector0
## Index vector indicating zero elements after rotation to df space
#  map.matrices=list(Full2DF=Full2DF,DF2full=DF2full,NormalVec2DF=NormalVec2DF,W=W)))
## 
###
#'  Spherization of multiway table array
#'
#' @param A An array or matrix of multiway table
#' @return tables A list of three tables: data :Input table A itself, expected : its expected table and diff: table of differecen between A and expected table
#' @return r.vec A integer vector indicating numbers of levels of A
#' @return df Degrees of freedom of A
#' @return matrices A list of four matrices, X,X.sub,X.inv and X.inv.sub, which are output of make.X.vector0, rotating to/from number-of-cell dimension and df dimension with/without meaningless rows/columns
#' @return zero.vec A vector indicating meaningful rows of rotations
#' @return map.matrices A list of four matrices: Full2df is a matrix transforming a table array vector to df-dimensional point, DF2full is a matrix transforming df-dimensional point to a table vector, Normalvec2DF is a matrix transforming test table to df-dimensional vector and W is a key matrix to calculate other matrices 
#' @examples
#' r <- c(2,3,4)
#' A <- array(1:prod(r),r)
#' t.sphere <- table.sphere(A)
table.sphere <- function(A){
  n <- sum(A) # sample size
  mg.cnt <- calc.marg(A) # marginal counts
  exp.table <- make.exp.table(A) # expected array
  diff.table <- make.diff.table(A) # original array - expected array
  chisq <- calc.chisq(A) # chisq value
  
  r.vec <- dim(A) # dimension vector	
  out <- make.X.vector0(r.vec) # zero element-indicating vector
  df <- sum(out$vector0) # degree of freedom
  W <- t(out$X.inv.sub) %*% diag(1/c(exp.table)) %*% (out$X.inv.sub)
  eigen.W <- eigen(W)
  tmp.diag <- tmp.diag.inv <- matrix(0,df,df)
  diag(tmp.diag) <- sqrt(eigen.W[[1]])
  diag(tmp.diag.inv) <- 1/sqrt(eigen.W[[1]])
  # Rotate an observe table vector to spherical space and its length^2 = chisq
  # |Full2DF %*% A|^2 = chisq
  Full2DF <- tmp.diag %*% t(eigen.W[[2]]) %*% out$X.sub
  # matrix rotates df-vector to full dimensional vector
  DF2full <- out$X.inv.sub %*% eigen.W[[2]] %*% tmp.diag.inv
  # Rotate test normal vector to spherical space normal vector
  NormalVec2DF <- out$X.inv.sub %*% (eigen.W[[2]]) %*% tmp.diag.inv  
  return(list(tables=list(data=A,expected=exp.table,diff=diff.table),r.vec = r.vec,df=df,matrices = list(X=out$X,X.sub=out$X.sub,X.inv=out$X.inv,X.inv.sub=out$X.inv.sub),zero.vec=out$vector0,map.matrices=list(Full2DF=Full2DF,DF2full=DF2full,NormalVec2DF=NormalVec2DF,W=W)))
}

###
# Arguments
## t.sphere: output of table.sphere()
## test.table: indicates alternative hypothesis direction in the shape of array
## odds.ratio.table: is an array that indicates the cells with which odds ratio should be set at the value or. (1,2,3,4) indicate cells as (cells[1] * cells[4])/(cells[2] * cells[3]) --> or
## n and K is arguments when serach the final result
# Values
## alt.table: Vectorized array
## alt.vec: Vector in df space
###
#'  Generate a table indicating alternative hypothesis
#'
#' @param t.sphere Output of table.sphere of a table of interest itself or a table sharing marginal counts with the table of interest
#' @param test.table A matrix or array of a table which indicates the direction of alternative hypothesis along with an argument odds.ratio.table
#' @param odds.ratio.table A matrix or array of a table with cell values being 0,1,2,3 or 4: Odds ratio of table is indicated by the (sum(cells of 1)*sum(cells of 4))/(sum(cells of 2)*sum(cells of 3)) 
#' @param or A number indicating odds ratio defined by odds.ratio.table argument
#' @param n Optional. An integer which makes a one-dimensional grid with n points, among which the table with odds ratio closest to the argument or is selected
#' @param K Optional. A number indicating the maximum value of square of chi-square value to search the table of alternative hypothsis
#' @return alt.table The estimated table of alternative hypothesis
#' @return alt.vec The vector of the estimated table in df-dimensional space
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t.sphere <- table.sphere(A)
#' test.table <- matrix(c(1,0.5,0,0,0,0),byrow=TRUE,2,3)
#' odds.ratio.table <- matrix(c(1,0,2,3,0,4),byrow=TRUE,2,3)
#' or <- 2
#' n <- 1000
#' K <- 10
#' make.alternative.table(t.sphere,test.table,odds.ratio.table,or,n,k)
make.alternative.table <- function(t.sphere,test.table,odds.ratio.table,or,n,K){
  # Default values
  if(missing(n))n<-1000
  if(missing(K))K<-10
  # alt.hy vector in df-space
  tmp <- c(test.table) %*% t.sphere$map.matrices$NormalVec2DF
  # dif vector in full-space
  dif <- t.sphere$map.matrices$DF2full %*% t(tmp)
  e.table <- t.sphere$tables$expected
  
  # Target OR should be calculated with XYZxyz
  X <- sum(e.table[which(odds.ratio.table==1)])
  Y <- sum(e.table[which(odds.ratio.table==2)])
  Z <- sum(e.table[which(odds.ratio.table==3)])
  W <- sum(e.table[which(odds.ratio.table==4)])
  x <- sum(dif[which(odds.ratio.table==1)])
  y <- sum(dif[which(odds.ratio.table==2)])
  z <- sum(dif[which(odds.ratio.table==3)])
  w <- sum(dif[which(odds.ratio.table==4)])
  # or of test.table
  tmp.or <- (X+x)*(W+w)/((Y+y)*(Z+z))
  # tmp.or and target or are compared and based on the comparison
  # search-range is defined
  if(tmp.or > or){
    ks <- seq(from=0,to=1,length=n)
  }else{
    ks <- seq(from=1,to=or/tmp.or*K,length=n)
  }
  # search-range is evenly spaced and their or is calculated
  # then closest one is selected.
  ors <- rep(0,n)
  for(i in 1:n){
    difs <- t.sphere$map.matrices$DF2full %*% t(tmp)*ks[i]
    
    x <- sum(difs[which(odds.ratio.table==1)])
    y <- sum(difs[which(odds.ratio.table==2)])
    z <- sum(difs[which(odds.ratio.table==3)])
    w <- sum(difs[which(odds.ratio.table==4)])
    ors[i] <- (X+x)*(W+w)/((Y+y)*(Z+z))
  }
  
  ret <- mean(ks[which(abs(ors-or)==min(abs(ors-or)))])
  ret2 <- t.sphere$map.matrices$DF2full %*% t(tmp)*ret
  return(list(alt.table=array(c(e.table)+ret2,t.sphere$r.vec),alt.vec=tmp*ret))
}
# tests is a list of tables representing test model and df-dimentional test vectors are generated using t.sphere information
#'  Generate test vectors from tables
#'
#' @param t.sphere table.sphere function's output
#' @param tests A list of arrays/matrix of tables representing tests
#' @return A number-of-test rows x df columns matrix of test vectos, each row of which is a test vector
#' @return alt.vec The vector of the estimated table in df-dimensional space
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t.sphere <- table.sphere(A)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' make.test.vecs(t.sphere,tests)

make.test.vecs <- function(t.sphere,tests){
  n.test <- length(tests)
  test.vecs <- matrix(0,n.test,t.sphere$df)
  for(i in 1:n.test){
    test.vecs[i,] <- c(c(tests[[i]]) %*% t.sphere$map.matrices$NormalVec2DF)
  }
  L.test.vecs <- sqrt(apply(test.vecs^2,1,sum))
  test.vecs/L.test.vecs
}
# df-dimensioanl unit random vectrors are generated and their inner product with test vectors (standardized statistic valus) area returned.
#'  Generate random K statistic values for test vectors under null hypothesis
#'
#' @param test.vecs A number-of-test rows x df columns matrix of test vectos, each row of which is a test vector
#' @param rmvn Optional A matrix of random unit vectors
#' @param n number of random unit vectors
#' @return A number-of-unit-vector rows x number-of-test columns matrix of K statistic values
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t.sphere <- table.sphere(A)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' test.vecs <- make.test.vecs(t.sphere,tests)
#' runitstat(test.vecs)
runitstat <- function(test.vecs,rmvn,n){
  # Default values
  if(missing(n))n <- 1000
  df <- length(test.vecs[1,])
  if(missing(rmvn)){
    rmvn <- st.mvn(n,df)
  }
  matrix((test.vecs %*% t(rmvn)),nrow=length(rmvn[,1]))
}
# When alternative hx is true, additional term on urunitstat is necessary to calculate statistic values. This function returns the coefficent for the modification.
#' Calculate K statistics of tests for output of make.alternative.table
#'
#' @param alt.tables Output of make.alternative.table, that is a list of table and vector in df-dimensional space of alternative hypothesis
#' @param test.vecs A number-of-test rows x df columns matrix of test vectos, each row of which is a test vector
#' @return A vector of K statisitc values
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t.sphere <- table.sphere(A)
#' test.table <- matrix(c(1,0.5,0,0,0,0),byrow=TRUE,2,3)
#' odds.ratio.table <- matrix(c(1,0,2,3,0,4),byrow=TRUE,2,3)
#' or <- 2
#' n <- 1000
#' K <- 10
#' alt.out <- make.alternative.table(t.sphere,test.table,odds.ratio.table,or,n,k)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' test.vecs <- make.test.vecs(t.sphere,tests)
#' alt.intersect(alt.out,test.vecs)
alt.intersect <- function(alt.tables,test.vecs){
  c(test.vecs %*% matrix(alt.tables$alt.vec,ncol=1))
}

# Of note, ks argument should be square-root of chisq values to be evaluated.
# Two terms are returned.
# term for null hx, depending on random vectors and test vectors. k/a,-k/a
# term for alternative hx, depending on test vectors. -b/a
#' Generate terms to calculate K statistic values of tests for random unit vectors
#'  
#' @param as NUmber-of-random-unit-vector rows x number-of-test matrix
#' @param ks A vector or K values (square of chi-square equivalent values) for which P-values should be estimated
#' @param bs Optional. A Vector of K values of alternative hypothesis table for tests
#' @return A list. k.per.a is 3-way array and is one set of terms to calculate K statistics, which is variable respect to unit vectors. b.per.a is another set of terms which K values for the alternative hypothesis itself
#' @examples
#' # utility function being used in other functions
window.ab <- function(as,ks,bs){
  zero.a <- which(as==0)
  as[zero.a] <- max(abs(as))*2
  # term for null hx, depending on random vectors and test vectors and cutoff test statistics. k/a,-k/a
  k.per.a <- outer(1/as,ks,"*")
  # term for alternative hx, depending on test vectors. -b/a
  if(missing(bs)){
    b.per.a <- NULL
  }else{
    #b.per.a <- -outer(1/as,bs,"*")
    b.per.a <- -t(t(1/as) * bs)
  }
  return(list(k.per.a=k.per.a,b.per.a=b.per.a))
}
# For each test, a segment in which |test statisitcs| is less than the cutoff stat value given after consideration of max stat among multiple tests.
# For each random vector, a segment is returned as a list of 1x2 matrix.
#' Calculate a segment to be integrated for individual K-cutoff value for each random unit vector
#'  
#' @param window.ab.out An output of window.ab
#' @param one.side Optional. Logical. When TRUE, test vectors indicate one-side test and otherwise two-sided
#' @return A list segments is a matrix with two rows, 1st column of which is the start and 2nd is the end of segments; ids is a vector indicating random unit vector ids to which segments belong 
#' @examples
#' # utility function being used in other functions
window.select <- function(window.ab.out,one.side){
  if(missing(one.side))one.side <- FALSE
  dm <- dim(window.ab.out$k.per.a)
  n.r <- dm[1]
  n.k <- dm[3]
  n.t <- dm[2]
  if(is.null(window.ab.out$b.per.a)){
    ret <- list()
    for(i in 1:n.k){
      tmp <- apply(abs(matrix(window.ab.out$k.per.a[,,i],ncol=n.t)),1,min)
      ret[[i]] <- list(segments=cbind(-tmp,tmp),ids=1:length(tmp))
    }
  }else{
    ret <- list()
    for(i in 1:n.k){
      ret[[i]] <- list()
      tmp.list.seg1 <- tmp.list.seg2 <- tmp.list.id <- list()
      for(j in 1:n.r){
        s <- cbind(c(-window.ab.out$k.per.a[j,,i]),c(window.ab.out$k.per.a[j,,i])) + window.ab.out$b.per.a[j,]
        if(one.side){
          s <- cbind(-c(window.ab.out$k.per.a[j,,i]),rep(0,length(window.ab.out$k.per.a[j,,i]))) + window.ab.out$b.per.a[j,]
        }
        tmp.seg <- segment.overlap(s)
        tmp.list.seg1[[j]] <- tmp.seg[1]
        tmp.list.seg2[[j]] <- tmp.seg[2]
        tmp.list.id[[j]] <- rep(j,length(tmp.list.seg1[[j]]))
      }
      ret[[i]] <- list(segments=cbind(unlist(tmp.list.seg1),unlist(tmp.list.seg2)),ids=unlist(tmp.list.id))
    }
  }
  ret
}
# for a list of segment, probability of df-chisq distribution within the segment is calculated and their mean is returned.
#' Calculate probability indicated by segments in random unit vector directions
#'  
#' @param segs.out Output of window.select
#' @param n NUmber of random unit vectors
#' @param df Degrees of freedom
#' @param one.side Optional. Logical. When TRUE, test vectors indicate one-side test and otherwise two-sided
#' @return Probability estimated
#' @examples
#' # utility function being used in other functions
calc.pr <- function(segs.out,n,df,one.side){
  if(missing(one.side))one.side <- FALSE
  ret <- rep(0,length(segs.out))
  for(i in 1:length(ret)){
    if(one.side){
      segs.out[[i]]$segments[which(segs<0)] <- 0
    }
    tmp1 <- pchisq.bid(segs.out[[i]]$segments[,1],df)
    tmp2 <- pchisq.bid(segs.out[[i]]$segments[,2],df)
    ret[i] <- (sum(tmp2)-sum(tmp1))/n
  }
  ret
}


# Power calculator
#' Power calculator
#'  
#' @param A An array or matrix of a table representing alternative hypothesis
#' @param tests A list of array/marix of tables indicating tests
#' @param alpha type 1 error threshold
#' @param n Number of iteration to estimate
#' @return 
#' \itemize{
#'  \item{"pvalues"}{null.p: estimated cumulative probability for null hypothesis when statistic K being x or higher 
#'  that is estimated for quantile of 1-alpha. alt.p: estimated value for alternative hypotheis }
#'  \item{"chisqs"}{Estimated K^2 value as a quantile of 1-alpha}
#'  \item{"alpha.given"}{Argument alpha}
#'  \item{"alpha.used"}{Estimated p.value for the estimated K^2 value}
#'  \item{"power"}{Power}
#' }
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' alpha <- 0.05
#' n <- 1000
#' power.mway.table(A,tests,alpha,n)
power.mway.table <- function(A,tests,alpha,n){
  if(missing(alpha))alpha <- 0.05
  if(missing(n))n<-1000
  tmp <- qmway.table(1-alpha,A,tests)
  x <- tmp$x
  alpha.estimated <- 1-tmp$p.estimated
  tmp.2 <- pmway.table(x,A,tests,nc=TRUE,n=n)
  null.p <- tmp.2$p.null
  alt.p <- tmp.2$p.alt
  #cut.off <- which(abs(null.p-(1-alpha))==min(abs(null.p-(1-alpha))))
  
  return(list(p.values = cbind(null.p,alt.p),chisqs = x,alpha.given=alpha,alpha.used=alpha.estimated,power=1-alt.p))
}
# With alpha-cut off, power is returned using two cumul prob for alt.hx.
#' Power calculator with statistic value being fixed
#'  
#' @param x K statistic values
#' @param rmway.out Output of rmway function
#' @param alpha Optional. type 1 error threshold
#' @return 
#' \itemize{
#'  \item{"pvalues"}{null.p: estimated cumulative probability for null hypothesis when statistic K being x or higher 
#'  that is estimated for quantile of 1-alpha. alt.p: estimated value for alternative hypotheis }
#'  \item{"chisqs"}{Estimated K^2 value as a quantile of 1-alpha}
#'  \item{"alpha.given"}{Argument alpha}
#'  \item{"alpha.used"}{Estimated p.value for the estimated K^2 value}
#'  \item{"power"}{Power}
#' }
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' n <- 1000
#' rmway.out <- rmway(n,A,tests)
#' x <- 8
#' alpha <- 0.05
#' power.mway(x,rmway.out,alpha)

power.mway <- function(x,rmway.out,alpha){
  if(missing(alpha))alpha <- 0.05
  null.p <- pmway(x,rmway.out,lower.tail=TRUE)
  alt.p <- pmway(x,rmway.out,alt.intersect=TRUE,lower.tail=TRUE)
  cut.off <- which(abs(null.p-(1-alpha))==min(abs(null.p-(1-alpha))))
  
  return(list(p.values = cbind(null.p,alt.p),chisqs = x^2,cut.off = x[cut.off]^2,power=1-alt.p[cut.off]))
}


# alt hx is given as intersect (b/a) instead of alt.hx table and cumul pr of chisq value x is returned.
#' Power calculator with intersect information rather than alternative hypothesis table given
#'  
#' @param x K statistic values
#' @param rmway.out Output of rmway function
#' @param intersect alternative hypothesis-dependent term
#' @param lower.tail ogical; if TRUE (default), probabilities are P[X ?????? x], otherwise, P[X > x].
#' @return Power
#' @examples
#' # Utility function for the authors
#' 
pmway.intersect <- function(x,rmway.out,intersect,lower.tail){
  if(missing(lower.tail))lower.tail <- TRUE
  w.ab.1 <- window.ab(rmway.out$runitstat.out,x,intersect)
  w.select.1 <- window.select(w.ab.1)
  ret <- calc.pr(w.select.1,rmway.out$n,rmway.out$t.sphere$df)
  if(!lower.tail){
    ret <- 1 - ret
  }
  ret
}

#' Generate althernative hypothesis tables in a direction with multipls K values
#'  
#' @param t.sphere Output of table.sphere
#' @param tests A list of tables indicating tests
#' @param test.table An array or matrix of table indicating the direction of alternaitve hypothesis
#' @param odds.ratio.table An array or matrix of table indicating 1,2,3,4 cells for calculation of odds ratio to be controled
#' @param ks Optional A vector of K statistic values
#' @return 
#' \itemize{
#' \item{"alt.tables"}{A list of alternative tables}
#' \item{"alt.vecs"}{A matrix of alternative table vectors in df-dimensional space}
#' \item{"alt.intersect.out"}{A matrix of terms used to calculate pertinent information}
#' \item{"stats"}{A vector of max K of alternative tables}
#' \item{"difs"}{A matrix each row of which has difference of cells of alternative hypothesis tables}
#' \item{"ors"}{A vector of odds ratios}
#' }
#' @examples
#' # Utility function for the authors
#'
model.table.series <- function(t.sphere,tests,test.table,odds.ratio.table,ks){
  if(missing(ks))ks <- seq(from=0,to=20,by=1)
  n.k <- length(ks)
  tmp <- c(test.table) %*% t.sphere$map.matrices$NormalVec2DF
  
  dif <- t.sphere$map.matrices$DF2full %*% t(tmp)
  e.table <- t.sphere$tables$expected
  test.vecs <- make.test.vecs(t.sphere,tests)
  difs <- matrix(0,n.k,length(e.table))
  alt.vecs <- matrix(0,n.k,t.sphere$df)
  stats <- ors <- rep(0,n.k)
  alt.tables <- list()
  alt.intersect.out <- matrix(0,n.k,length(test.vecs[,1]))
  
  X <- sum(e.table[which(odds.ratio.table==1)])
  Y <- sum(e.table[which(odds.ratio.table==2)])
  Z <- sum(e.table[which(odds.ratio.table==3)])
  W <- sum(e.table[which(odds.ratio.table==4)])
  
  for(i in 1:n.k){
    alt.vecs[i,] <- tmp * ks[i]
    alt.intersect.out[i,] <- c(test.vecs %*% matrix(alt.vecs[i,],ncol=1))
    stats[i] <- max((test.vecs %*% matrix(alt.vecs[i,],ncol=1))^2)
    alt.tables[[i]] <- array(c(e.table)+ t.sphere$map.matrices$DF2full %*% matrix(alt.vecs[i,],ncol=1),t.sphere$r.vec)
    difs[i,] <- t.sphere$map.matrices$DF2full %*% matrix(alt.vecs[i,],ncol=1)
    tmp.dif <- difs[i,]
    x <- sum(tmp.dif[which(odds.ratio.table==1)])
    y <- sum(tmp.dif[which(odds.ratio.table==2)])
    z <- sum(tmp.dif[which(odds.ratio.table==3)])
    w <- sum(tmp.dif[which(odds.ratio.table==4)])
    ors[i] <- (X+x)*(W+w)/((Y+y)*(Z+z))
  }
  
  return(list(alt.tables=alt.tables,alt.vecs=alt.vecs,alt.intersect.out=alt.intersect.out,stats=stats,difs=difs,ors=ors))
  
}
#' Generate random tables under unll or alternative hypothesis
#'  
#' @param n Number of random tables to be generated
#' @param A An array or matrix of table that specifies marginal counts and also alternative hypothesis if nc is TRUE
#' @param tests Optional. When TRUE, test statistics and maximum among them are calculated and returned
#' @param nc Optional. When FALSE, random tables are in null hypothesis and when TRUE in alternative hypothesis with its ccenter being table A
#' @param one.side Optional. When TRUE, tests are one-sided, otherwise two-sided 
#' @return 
#' \itemize{
#' \item{"stats.all"}{A matrix of K statistics of random tables}
#' \item{"stats"}{A matrix of K^2 statistics of random tables}
#' \item{"rtables"}{A list of random tables}
#' \item{"ellipsoid.loc"}{A matrix of coordinates of random tables that are in df-dimensional space but the contours are ellipsoids}
#' }
#' @details This function is one of three major functions in this package. Functions for a distribution are rxxxx,dxxxx, pxxxx and qxxxx
#' in general and this package provides three of them, rmway.table, pmway.table and qmway.table functions.
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' n <- 1000
#' rmway.tables.null <- rmway.table(n,A)
#' rmway.tables.alt <- rmway.table(n,A,nc=TRUE)
#'
rmway.table <- function(n,A,tests,nc,one.side){
  t.sphere <- table.sphere(A)
  if(missing(tests))tests <- NULL
  if(missing(nc))nc <- FALSE
  if(missing(one.side))one.side <- FALSE
  
  if(!nc){
    alt.vec <- rep(0,t.sphere$df)
  }else{
    alt.vec <- t.sphere$map.matrices$Full2DF %*% c(t.sphere$tables$diff)
  }
  if(!is.null(tests)){
    test.vecs <- make.test.vecs(t.sphere,tests)
  }
  
  #rmvn <- st.mvn(n,t.sphere$df)
  rmvn <- matrix(rnorm(n*t.sphere$df),nrow=n)
  loc.df <- t(rmvn) + c(alt.vec)
  #loc.full <- t.sphere$matrices$X.inv.sub %*% loc.df
  loc.full <- t.sphere$map.matrices$DF2full %*% loc.df
  ellipsoid.loc <- t.sphere$matrices$X.sub %*% loc.full
  loc.full.2 <- loc.full + c(t.sphere$tables$expected)
  rtables <- list()
  for(i in 1:n){
    rtables[[i]] <- array(loc.full.2[,i],t.sphere$r.vec)
  }
  if(is.null(tests)){
    stats <- NULL
    stats.all <- NULL
  }else{
    stats <- test.vecs %*% loc.df
    if(one.side){
      stats[which(stats<0)] <- 0
    }
    stats.all <- stats^2
    stats <- apply(stats,2,max)
  }
  
  return(list(stats.all = stats.all,stats=stats^2,rtables=rtables,ellipsoid.loc=t(ellipsoid.loc)))
}
#' Estimate p-values/cumulative probability of maximum statistics on multi-way table under unll or alternative hypothesis
#'  
#' @param x A vector of K^2 quantiles
#' @param A An array or matrix of table that specifies marginal counts
#' @param tests A list of tables
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ?????? x], otherwise, P[X > x].
#' @param nc Optional. if FALSE (default), random tables are in null hypothesis and otherwise in alternative hypothesis with its ccenter being table A
#' @param one.side Optional. When TRUE, tests are one-sided, otherwise two-sided 
#' @param n number of random unit vectors used for estimation
#' @return 
#' \itemize{
#' \item{"p.null"}{cumulative probability under null hypothesis}
#' \item{"p.alt"}{cumulative probability under alternative hypothesis}
#' }
#' @details This function is one of three major functions in this package. Functions for a distribution are rxxxx,dxxxx, pxxxx and qxxxx
#' in general and this package provides three of them, rmway.table, pmway.table and qmway.table functions.
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' n <- 1000
#' x <- 0:10
#' pmway.table(x,A,tests)
#'
pmway.table <- function(x,A,tests,lower.tail,nc,one.side,n){
  if(missing(lower.tail))lower.tail <- TRUE
  if(missing(nc))nc <- FALSE
  if(missing(one.side))one.side <- FALSE
  if(missing(n))n <- 1000
  if(!nc){
    ret.null <- pmway.table.null(x,A,tests,lower.tail=lower.tail,one.side=one.side,n=n)
    ret.alt <- NULL
    return(list(p.null=ret.null,p.alt=ret.alt))
  }
  t.sphere <- table.sphere(A)
  if(!nc){
    alt.vec <- rep(0,t.sphere$df)
  }else{
    
    alt.vec <- t.sphere$matrices$X.sub %*% c(A)
  }
  test.vecs <- make.test.vecs(t.sphere,tests)
  rmvn <- st.mvn(n,t.sphere$df)
  runitstat.out <- runitstat(test.vecs,rmvn)
  intersect <- c(test.vecs %*% matrix(alt.vec,ncol=1))
  w.ab.1 <- window.ab(runitstat.out,sqrt(x),intersect)
  w.select.1 <- window.select(w.ab.1,one.side=one.side)
  pr.sum <- calc.pr(w.select.1,n,t.sphere$df)
  
  max.stat <- apply(abs(runitstat.out),1,max)
  non.zero <- which(max.stat !=0)
  k.per.a <- outer(1/max.stat,sqrt(x),"*")
  pr.null <- pchisq(k.per.a^2,t.sphere$df)
  pr.sum.null <- apply(pr.null,2,sum)/n
  
  if(lower.tail){
    if(one.side){
      ret <- 1 - (1 - pr.sum)/2
      ret.null <- 1 - (1 - pr.sum.null)/2
    }else{
      ret <- pr.sum
      ret.null <- pr.sum.null
    }
  }else{
    if(one.side){
      ret <- (1 - pr.sum)/2
      ret.null <- (1 - pr.sum.null)/2
    }else{
      ret <- 1 - pr.sum
      ret.null <- 1 - pr.sum.null
    }
  }
  return(list(p.null=ret.null,p.alt=ret))
}


# This is important.
#' Estimate p-values/cumulative probability of maximum statistics on multi-way table under unll hypothesis
#'  
#' @param x A vector of K^2 quantiles
#' @param A An array or matrix of table that specifies marginal counts
#' @param tests A list of tables
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ?????? x], otherwise, P[X > x].
#' @param one.side Optional. When TRUE, tests are one-sided, otherwise two-sided 
#' @param n number of random unit vectors used for estimation
#' @return cumulative probability under null hypothesis
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' n <- 1000
#' x <- 0:10
#' pmway.table.null(x,A,tests)
pmway.table.null <- function(x,A,tests,lower.tail,one.side,n){
  if(missing(lower.tail))lower.tail <- TRUE
  if(missing(one.side))one.side <- FALSE
  if(missing(n))n <- 1000
  t.sphere <- table.sphere(A)
  test.vecs <- make.test.vecs(t.sphere,tests)
  rmvn <- st.mvn(n,t.sphere$df)
  runitstat.out <- runitstat(test.vecs,rmvn)
  max.stat <- apply(abs(runitstat.out),1,max)
  #non.zero <- which(max.stat !=0)
  k.per.a <- outer(1/max.stat,sqrt(x),"*")
  pr <- pchisq(k.per.a^2,t.sphere$df)
  pr.sum <- apply(pr,2,sum)/n
  
  if(lower.tail){
    if(one.side){
      ret <- 1 - (1 - pr.sum)/2
    }else{
      ret <- pr.sum
    }
  }else{
    if(one.side){
      ret <- (1 - pr.sum)/2
    }else{
      ret <- 1 - pr.sum
    }
  }
  ret
}

# This is important to construct qmway.table()
#' Grossly estimate quantile value. Utility function for qmway.table function
#'  
#' @param p probability
#' @param A An array or matrix of table that specifies marginal counts
#' @param tests A list of tables
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ?????? x], otherwise, P[X > x].
#' @param x a vector of K statistic values
#' @param nc if FALSE (default) null hypothesis, otherwise alternative hypothesis
#' @param one.side Optional. When TRUE, tests are one-sided, otherwise two-sided 
#' @param n number of random unit vectors used for estimation
#' @return two quantile values and two corresponding p values between which p and its quantile exist
#' @examples
#' # Utility function
qmway.table.pre <- function(p,A,tests,lower.tail,x,nc,one.side,n){
  if(missing(lower.tail))lower.tail <- TRUE
  if(missing(x))x <- seq(from=0,to=100,by=10)
  if(missing(nc))nc <- FALSE
  if(missing(one.side))one.side <- FALSE
  if(missing(n))n <- 1000
  if(!lower.tail){
    p <- 1-p
  }
  p.out <- pmway.table(x,A,tests,lower.tail=lower.tail,nc=nc,one.side=one.side,n=n)$p.null
  if(p-min(p.out) < 0){
    print("min of x is too big")
    return(NULL)
  }else if(max(p.out) - p < 0){
    print("max of x is too small")
    return(NULL)
  }
  tmp.1 <- which(p.out <= p)
  tmp.1 <- tmp.1[length(tmp.1)]
  tmp.2 <- which(p.out >= p)
  tmp.2 <- tmp.2[1]
  
  return(c(x[tmp.1],x[tmp.2],p.out[tmp.1],p.out[tmp.2]))
}
# qmway.table returns stat value for p given as an argument. The returned value is selected by interporating between the closest stat values by two step interporations.
#' Estimate quantiles of maximum statistics on multi-way table under unll or alternative hypothesis
#'  
#' @param p probability vector
#' @param A An array or matrix of table that specifies marginal counts
#' @param tests A list of tables
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ?????? x], otherwise, P[X > x].
#' @param x a vector of K statistic values
#' @param nc Optional. if FALSE (default), random tables are in null hypothesis and otherwise in alternative hypothesis with its ccenter being table A
#' @param one.side Optional. When TRUE, tests are one-sided, otherwise two-sided 
#' @param n number of random unit vectors used for estimation
#' @return 
#' \itemize{
#' \item{"x"}{quantile}
#' \item{"p.used"}{Argument p}
#' \item{"p.estimated"}{p value corresponding to the estimated quantile, which is some deviated from p of argument because, estimation process is involeved.}
#' }
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' n <- 1000
#' p <- 0.05
#' qmway.table(p,A,tests)
#'
qmway.table <- function(p,A,tests,lower.tail,x,nc,one.side,n){
  if(missing(lower.tail))lower.tail <- TRUE
  if(missing(x))x <- seq(from=0,to=100,by=10)
  if(missing(nc))nc <- FALSE
  if(missing(one.side))one.side <- FALSE
  if(missing(n))n <- 1000
  pre.out <- qmway.table.pre(p,A,tests,lower.tail=lower.tail,x=x,nc=nc,one.side=one.side,n=n)
  post.x <- seq(from=pre.out[1],to=pre.out[2],length=n)
  post.out <- qmway.table.pre(p,A,tests,lower.tail=lower.tail,x=post.x,nc=nc,one.side=one.side,n=n)
  ret <- post.out[1] + (post.out[2]-post.out[1]) * (p-post.out[3])/(post.out[4]-post.out[3])
  ret.p <- pmway.table(ret,A,tests)$p.null
  return(list(x=ret,p.used=p,p.estimated=ret.p))
}

# Tests
#' Max test 
#'  
#' @param A An array or matrix of table that specifies marginal counts
#' @param tests A list of tables
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ?????? x], otherwise, P[X > x].
#' @param one.side Optional. When TRUE, tests are one-sided, otherwise two-sided 
#' @param n number of random unit vectors used for estimation
#' @return 
#' \itemize{
#' \item{"statistic"}{maximum among test statistics; K^2-equivalent}
#' \item{"p.value"}{multiple testing-corrected P value}
#' }
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' n <- 1000
#' mway.test(A,tests,n=10000)
mway.test <- function(A,tests,lower.tail,one.side,n){
  if(missing(lower.tail))lower.tail <- FALSE
  if(missing(one.side))one.side <- one.side <- FALSE
  if(missing(n))n <- 1000
  t.sphere <- table.sphere(A)
  test.vecs <- make.test.vecs(t.sphere,tests)
  stats <- test.vecs %*% t.sphere$map.matrices$Full2DF %*% c(t.sphere$tables$diff)
  
  max.stats <- max(abs(stats))
  if(one.side){
    max.stats <- max(stats)
  }
  p <- pmway.table(max.stats^2,A,tests,lower.tail=lower.tail,one.side=one.side,n=n)
  return(list(statistic=max.stats^2,p.value=p))
}
#' Generate data matrix from table
#' 
#' @param m An array or matrix of table that specifies marginal counts
#' @return A matrix of sample size x number of variables, levels of which are integers starting from 1
#' @examples
#' m <- matrix(c(10,20,30,40,50,60),2,3)
#' make.data.from.table(m)
make.data.from.table <- function(m){
  r <- dim(m)
  address <- which(m <= Inf,arr.ind=TRUE)
  N <- sum(m)
  x <- matrix(0,0,length(r))
  for(i in 1:length(address[,1])){
    tmp <- matrix(rep(address[i,],m[i]),byrow=TRUE,nrow=m[i])
    x <- rbind(x,tmp)
  }
  x
}
#' Generate an array or matrix from data matrix
#'  
#' @param d a matrix of sample size x number of variables, levels of which should be integers starting from 1
#' @return An array or matrix
#' @examples
#' m <- matrix(c(10,20,30,40,50,60),2,3)
#' d <- make.data.from.table(m)
#' make.table.from.data(d)
make.table.from.data <- function(d){
  rg <- apply(d,2,max)
  ret <- array(rep(0,prod(rg)),rg)
  prod.rg <- cumprod(rg)
  prod.rg <- c(1,prod.rg[-length(rg)])
  tmp <- apply(t(d-1) * prod.rg,2,sum) + 1
  for(i in 1:length(tmp)){
    ret[tmp[i]] <- ret[tmp[i]]+1
  }
  ret
}
#' Randomize a data matrix
#'  
#' @param data a matrix of sample size x number of variables
#' @param v a vector of 0 or 1 with length being the number of variables, which indicates variables to be shuffled.
#' @return data matrix 
#' @examples
#' m <- matrix(c(10,20,30,40,50,60),2,3)
#' d <- make.data.from.table(m)
#' make.table.from.data(d)
perm.data <- function(data,v){
  if(missing(v))v <- rep(1,length(data[1,]))
  ret <- data
  for(i in 1:length(data[1,])){
    if(v[i]==1){
      ret[,i] <- sample(ret[,i])
    }
  }
  ret
}



# rmway.out has info of alt hx table and otherts. 
# x is chisq value and p values (cumulative prob) is returned
# When alt.intersect is TRUE, cumulative prob of x is returned, otherwise cumulative prob of x when null is true.
#' Randomize a data matrix
#'  
#' @param x K statistic value
#' @param rmway.out output of rmway function
#' @param alt.intersect logical, if false (default), null hypotheis is assumued, otherwise alternative hypothesis is assumed
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ?????? x], otherwise, P[X > x].
#' @return probability
#' @examples
#' # Utility function
pmway <- function(x,rmway.out,alt.intersect,lower.tail){
  if(missing(alt.intersect))alt.intersect <- FALSE
  if(missing(lower.tail))lower.tail <- TRUE
  if(!alt.intersect){
    w.ab.1 <- window.ab(rmway.out$runitstat.out,x)
  }else if(alt.intersect){
    w.ab.1 <- window.ab(rmway.out$runitstat.out,x,rmway.out$alt.intersect.out)
  }
  w.select.1 <- window.select(w.ab.1)
  ret <- calc.pr(w.select.1,rmway.out$n,rmway.out$t.sphere$df)
  if(!lower.tail){
    ret <- 1 - ret
  }
  ret
}
# rmway.out is already calculated and given as an argment and alt.hx is given to calculate alt hx table with other info.
#' Generate random K statistic values with rmway function's output reused
#'  
#' @param alt.hx A matrix or array of a table which indicates the direction of alternative hypothesis along with an argument odds.ratio.table
#' @param odds.ratio.table A matrix or array of a table with cell values being 0,1,2,3 or 4: Odds ratio of table is indicated by the (sum(cells of 1)*sum(cells of 4))/(sum(cells of 2)*sum(cells of 3)) 
#' @param or A number indicating odds ratio defined by odds.ratio.table argument
#' @param rmway.out rmway function's output
#' @param n.alt Optional. An integer which makes a one-dimensional grid with n points, among which the table with odds ratio closest to the argument or is selected
#' @param K.alt Optional. A number indicating the maximum value of square of chi-square value to search the table of alternative hypothsis
#' @return 
#' \itemize{
#'  \item{"t.sphere"}{ Output of table.sphere(A)}
#'  \item{"alt.table"}{Argument alt.table itself}
#'  \item{"tests"}{Argument tests}
#'  \item{"test.vecs"}{Matrix whose rows are test vectors}
#'  \item{"runitstat.out"}{Output of runitstat function}
#'  \item{"alt.intersect.out"}{Output of alt.intersect function}
#'  \item{"n"}{Argument n}
#'  \item{"A"}{Argument A}
#'  \item{"alt.hx"}{Argument alt.hx}
#'  \item{"odds.ratio.table"}{Argument odds.ratio.table}
#'  \item{"or"}{Argument or}
#'  \item{"n.alt"}{Argument n.alt}
#'  \item{"K.alt"}{Argument K.alt}
#'  \item{"rmvn"}{A matrix of random unit vectors, each row of which is a random unit vector}
#'  }
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' n <- 1000
#' rmway.out <- rmway(n,A,tests)
#' test.table <- matrix(c(1,0.5,0,0,0,0),byrow=TRUE,2,3)
#' odds.ratio.table <- matrix(c(1,0,2,3,0,4),byrow=TRUE,2,3)
#' or <- 2
#' rmway.alt.only(test.table,odds.ratio.table,or,rmway.out)
rmway.alt.only <- function(alt.hx,odds.ratio.table,or,rmway.out,n.alt,K.alt){
  if(missing(n.alt))n.alt <- 1000
  if(missing(K.alt))K.alt <- 10
  t.sphere <- rmway.out$t.sphere
  alt.table <- make.alternative.table(t.sphere,alt.hx,odds.ratio.table,or,n.alt,K.alt)
  test.vecs <- rmway.out$test.vecs
  rmvn <- rmway.out$rmvn
  runitstat.out <- rmway.out$runitstat.out
  alt.intersect.out <- alt.intersect(alt.table,test.vecs)
  return(list(t.sphere=t.sphere,alt.table=alt.table,tests=tests,test.vecs=test.vecs,runitstat.out=runitstat.out,alt.intersect.out=alt.intersect.out,n=n,A=rmway.out$tables$data,tests=tests,alt.hx=alt.hx,odds.ratio.table=odds.ratio.table,or=or,n.alt=n.alt,K.alt=K.alt,rmvn = rmvn))
}
# alt hx is defined with alt.hx being an array and target odds ratio defined by odds ratio table with {1,2,3,4} with target or value. n.alt and K.alt are parameters to select alt hx table.
# alt hx table is returned with other infos.
#' Generate random K statistic values
#'  
#' @param n Number of random unit vectors
#' @param A An array or matrix of a table
#' @param tests A list of test tables
#' @param alt.hx A matrix or array of a table which indicates the direction of alternative hypothesis along with an argument odds.ratio.table
#' @param odds.ratio.table A matrix or array of a table with cell values being 0,1,2,3 or 4: Odds ratio of table is indicated by the (sum(cells of 1)*sum(cells of 4))/(sum(cells of 2)*sum(cells of 3)) 
#' @param or A number indicating odds ratio defined by odds.ratio.table argument
#' @param n.alt Optional. An integer which makes a one-dimensional grid with n points, among which the table with odds ratio closest to the argument or is selected
#' @param K.alt Optional. A number indicating the maximum value of square of chi-square value to search the table of alternative hypothsis
#' @return 
#' \itemize{
#'  \item{"t.sphere"}{ Output of table.sphere(A)}
#'  \item{"alt.table"}{Argument alt.table itself}
#'  \item{"tests"}{Argument tests}
#'  \item{"test.vecs"}{Matrix whose rows are test vectors}
#'  \item{"runitstat.out"}{Output of runitstat function}
#'  \item{"alt.intersect.out"}{Output of alt.intersect function}
#'  \item{"n"}{Argument n}
#'  \item{"A"}{Argument A}
#'  \item{"alt.hx"}{Argument alt.hx}
#'  \item{"odds.ratio.table"}{Argument odds.ratio.table}
#'  \item{"or"}{Argument or}
#'  \item{"n.alt"}{Argument n.alt}
#'  \item{"K.alt"}{Argument K.alt}
#'  \item{"rmvn"}{A matrix of random unit vectors, each row of which is a random unit vector}
#'  }
#' @examples
#' A <- matrix(c(10,20,30,40,50,60),2,3)
#' t1 <- t2 <- t3 <- matrix(c(1,0,0,0,0,0),byrow=TRUE,2,3)
#' t2[1,2] <- 0.5
#' t3[1,2] <- 1
#' tests <- list(t1,t2,t3)
#' n <- 1000
#' rmway.out <- rmway(n,A,tests)
rmway <- function(n,A,tests,alt.hx,odds.ratio.table,or,n.alt,K.alt){
  
  t.sphere <- table.sphere(A)
  if(missing(alt.hx)){
    alt.hx <- odds.ratio.table <- or <- n.alt <- K.alt <- NULL
    alt.table <- list(alt.table=t.sphere$tables$expected,alt.vec=rep(0,t.sphere$df))
  }else{
    if(missing(n.alt))n.alt <- 1000
    if(missing(K.alt))K.alt <- 10
    alt.table <- make.alternative.table(t.sphere,alt.hx,odds.ratio.table,or,n.alt,K.alt)
  }
  test.vecs <- make.test.vecs(t.sphere,tests)
  rmvn <- st.mvn(n,t.sphere$df)
  runitstat.out <- runitstat(test.vecs,rmvn)
  alt.intersect.out <- alt.intersect(alt.table,test.vecs)
  return(list(t.sphere=t.sphere,alt.table=alt.table,tests=tests,test.vecs=test.vecs,runitstat.out=runitstat.out,alt.intersect.out=alt.intersect.out,n=n,A=A,tests=tests,alt.hx=alt.hx,odds.ratio.table=odds.ratio.table,or=or,n.alt=n.alt,K.alt=K.alt,rmvn = rmvn))
}
# A??????????????????0??????????????????????????????????????????(???????????????)
# A???????????????????????????(?????????????????????)
# A.ori???NULL???????????????????????????A.ori??????????????????????????????????????????????????????A??????????????????????????????
#' Utility function to handle tables with cells whose expected value is 0
#'
#' @param A an array or matrix of table, that may be the table with cells whose expected value is zero (0 cells) or 
#' that may be the table without such cells
#' @param A.ori if missing(default), A represetnts a table with 0 cell(s) and the function converts A to a table without 0 cell(s).
#' Otherwise, A should be a table without 0 cell(s) and A.ori is the original array/matrix before shrinkage.
#' of table with 0 cells
#' @return  An array or matrix of table. if A.ori is missing, shrunken table is returned otherwise a table with 0 cells recovered.
#' with 0 cells.
#' @examples
#' n <- 5
#' M <- list()
#' N <-100
#' for(i in 1:n){
#'  M[[i]] <- sample(0:5,sample(2:4,1))
#'  M[[i]] <- M[[i]]/sum(M[[i]])*N
#' }
#' tmp <- M[[1]]
#' for(i in 2:length(M)){
#'   tmp <- t(M[[i]]/N) %x% tmp
#' }
#' A <- array(tmp,lapply(M,length))
#'
#' A.r <- reshape.table(A)
#' A.2 <- reshape.table(A.r,A)
#' calc.marg(A)
#' calc.marg(A.r)
#' calc.marg(A.2)
#' A-A.2 
#' 
reshape.table <- function(A,A.ori){
  if(missing(A.ori)){
    A.ori <- NULL
  }
	inv <- FALSE
	if(!is.null(A.ori)){
		inv <- TRUE
	}
	if(!inv){
		A.ori <- A
	}
	M <- calc.marg(A.ori)
	S <- lapply(M,sign)
	new.dim <- lapply(S,sum)
	tmp <- S[[1]]
	for(i in 2:length(S)){
		tmp <- t(S[[i]]) %x% tmp
	}
	if(!inv){
		return(array(A[which(tmp!=0)],new.dim))
	}else{
		tmp.A <- rep(0,length(A.ori))
		tmp.A[which(tmp!=0)] <- c(A)
		return(array(tmp.A,dim(A.ori)))
	}	
}
