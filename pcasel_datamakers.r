library(MASS)
library(mvtnorm)


#' Data generator for genotypes based on mvn
#'
#' Here our aim is to generate random genotypes by
#' first simualting an allele frequnecy p using a
#' multivariate normal, and then sampling a binomial
#' Binom(2,p) genotype. I simulate a population by first
#' generating a covmatrix from points in the square
#' [0:1] x [0:1], then using 2 - dist(x,y) as cov
#' matrix. To simulate selection, I move the points in 
#' [0:sr,0:sr] to [se:sr+se;se:sr+se], where sr and se
#' are the selection region and selection strength,
#' respectively, simulating samples "further" away.
#' 
#' @param n.samples the number of sampled individuals
#' @param n.neutral.snps number of neutral snps
#' @param n.selected.snps number of snps with selection
datamaker.mvngenotypes <- function(n.samples=100, 
                                   n.neutral.snps=9900,
                                   n.selected.snps=0100,
                                   sel.region=sqrt(0.25),
                                   sel.strength=4,...){

    pts.x <- runif(n.samples)
    pts.y <- runif(n.samples)
    pts <- cbind(pts.x,pts.y)
    pts <- pts[order(pts.x,pts.y),]
    dmat <- dist(pts)
    Sigma.neutral <- max(dmat) - as.matrix(dmat)
    
    which.sel <- pts.x < sel.region & pts.y < sel.region
    pts[which.sel,] <- pts[which.sel,] - sel.strength
    dmat <- dist(pts)
    Sigma.sel <- max(dmat) - as.matrix(dmat)
    
    mu <- rep(0.5, n.samples)
    
    data.neutral <- mvrnorm(n.neutral.snps,
                            mu=mu,
                            Sigma=Sigma.neutral)
    data.sel <- mvrnorm(n.selected.snps,
                            mu=mu,
                            Sigma=Sigma.sel)
    print("simulated data")
    
    data <- rbind(data.neutral, data.sel)


    meta <- list()
    meta$selpop <- which.sel
    meta$issel <- c(rep(F, n.neutral.snps),
                    rep(T, n.selected.snps))
    
    return(list(input=data, meta=meta))
}

#' dataset simulator
#' 
#' this is the main function to generate data sets.
#' It calls different other function to do specific 
#' simulations.
#' @param args A list of arguments
#' @param return a list 
datamaker <- function(args){
    data <- do.call(args$fun, args)
    if(!is.null(args$geno) && args$geno){
        spread <- 1
        if(!is.null(args$spread)){
            spread <- args$spread
        }
        print(spread)
        if(length(spread)==2){
            data$input <- mvn2snp.random(data$input, spread)
        }else{
	data$input <- mvn2snp.basic(data$input, spread)
        }
    }
    raw <- data$input
    data$input <- list()
    data$input$raw <- raw
    data$input$svd <- svd(centerscale(raw))
    return(data)
}

mvn2snp.random <- function(input, spread=c(0.05,0.5), 
                           mean.freq.mean=0.3, mean.freq.sd=0.1){
    median.sd <- median(apply(input, 1, sd))
    n <- nrow(input)
    rspread <- runif(n,spread[1],spread[2])
    mean.freq <- rnorm(n,mean.freq.mean, mean.freq.sd)

    i2 <- (input/median.sd)*spread +mean.freq
    i3 <- pmin(pmax(i2,0),1)
    i4 <- matrix(rbinom(prod(dim(input)), 2, i3),
		 ncol=ncol(input))
    return(i4)
}


mvn2snp.basic <- function(input, spread=2){
    median.sd <- median(apply(input, 1, sd))
    mean.freq <- 0.5

    i2 <- (input/median.sd)*spread +mean.freq
    i3 <- pmin(pmax(i2,0),1)
    i4 <- matrix(rbinom(prod(dim(input)), 2, i3),
		 ncol=ncol(input))
    return(i4)
}

datamaker.simulate.from.basis <- function(n.snps, basis, evs, 
                                          Sigma.to.add=0, ...){
    #evs <- c(1,1, rep(0,length(evs)-2))
    Sigma <- basis %*% diag(evs) %*% t(basis)
    print(min(evs))
    Sigma <- Sigma + Sigma.to.add

    input <- rmvnorm(n.snps, sigma=Sigma, 
		     method="svd")
}

datamaker.discrete.cosine.peaksel <- function(
				      n.neutral.snps=9500,
				      n.selected.snps=0500,
				      n.samples.per.dim=10,
				      neutral.ev.coef=c(2,1),
				      sel.pos = c(0.5,0.5),
                                      sel.scale = 4,
                                      sel.var = 0.1,
                                      sel.cov = 0,
				      ...){
    nspd = n.samples.per.dim
    g <- discrete.cosine.basis.2d(nspd, nspd, nspd,
				  t=neutral.ev.coef, sort=T)
    g$evs <- pmax(g$evs,1e-6)
    selbasis <- datamaker.selbasis(n.samples.per.dim=n.samples.per.dim,
                                   sel.pos, sel.var, sel.cov,
                                   sel.scale)
    #selbasis <- selbasis[g$order]
    input <- datamaker.simulate.from.basis(n.neutral.snps,
                                           g$basis, g$evs)
    input.sel <- datamaker.simulate.from.basis(n.selected.snps,
                                           g$basis, g$evs, selbasis)
    meta <- list()
    meta$issel <- c(rep(F, n.neutral.snps),
                    rep(T, n.selected.snps))

    return(list(input=rbind(input, input.sel),
		meta=meta))
}

datamaker.discrete.cosine <- function(
				      n.neutral.snps=9500,
				      n.selected.snps=0500,
				      n.samples.per.dim=10,
				      neutral.ev.coef=c(2,1),
				      sel.ev.coef=c(3,1),
				      ...){
    nspd = n.samples.per.dim
    g <- discrete.cosine.basis.2d(nspd, nspd, nspd,
				  t=neutral.ev.coef)
    g$evs <- pmax(g$evs,1e-6)
    input <- datamaker.simulate.from.basis(n.neutral.snps,
                                           g$basis, g$evs)

    evs.sel <- g$evs
    nsel <- length(sel.ev.coef)
    sel.ev.coef <- 1 + sel.ev.coef #additive selection
    evs.sel[2:(nsel+1)] <- evs.sel[2:(nsel+1)]*sel.ev.coef

    input.sel <- datamaker.simulate.from.basis(n.selected.snps,
                                               g$basis, evs.sel)

    meta <- list()
    meta$issel <- c(rep(F, n.neutral.snps),
                    rep(T, n.selected.snps))

    return(list(input=rbind(input, input.sel),
		meta=meta))
}

datamaker.discrete.cosine2 <- function(x.scale, y.scale,
                                       sel.strength,
                                       sel.angle,
                                       dry=F,
                                       ...){
    neutral.ev.coef <- c(x.scale, y.scale)
    sel.x <- sin(sel.angle) * sel.strength * x.scale
    sel.y <- cos(sel.angle) * sel.strength * y.scale
    sel.ev.coef <- c(sel.x, sel.y)

    if(dry) return(cbind(neutral.ev.coef, sel.ev.coef))
    datamaker.discrete.cosine(..., neutral.ev.coef=neutral.ev.coef,
                              sel.ev.coef=sel.ev.coef
                              )
}

datamaker.selbasis <- function(n.samples.per.dim=10,
                                                mvn.mean=c(.5,.5),
                                                mvn.var=c(.5,.5),
                                                mvn.cov=0,
                                                scale=4
                                                ){
    v <- (1:n.samples.per.dim-1)/(n.samples.per.dim-1)
    g <- expand.grid(v,v) 
    Sigma <- matrix(NA, ncol=2,nrow=2)
    diag(Sigma) <- mvn.var
    Sigma[1,2] <- mvn.cov
    Sigma[2,1] <- mvn.cov

    mvd <- dmvnorm(g, mvn.mean, Sigma)
    mvd <- mvd %*% t(mvd)
    mvd <- mvd - min(mvd)
    mvd <- mvd / max(mvd)
    return(mvd * scale)
}



discrete.cosine.basis.1d <- function(n=50,n.vectors=n){
    sapply(0:n.vectors,function(k)cos((0:n+0.5)*k*pi/n))
}
discrete.cosine.basis.2d <- function(n=20,m=20, n.vectors=n,
				     t=c(1,1), sort=T){
    grid <- expand.grid(0:(n.vectors-1),0:(n.vectors-1))
    if(sort){
        o <- order(rowSums(grid),grid[,1])
        grid <- grid[o,]
    }

    basis <- apply(grid,1, function(k)
		 discrete.cosine.2d(n,m,k))
    evs <- discrete.cosine.evs(n.vectors, t)

    if(sort)
        return(list(basis=basis, evs=evs, grid=grid, order=o))
    return(list(basis=basis, evs=evs, grid=grid))
    
}

discrete.cosine.evs <- function(n.vectors,t){
    grid <- expand.grid(0:(n.vectors-1),0:(n.vectors-1))
    o <- order(rowSums(grid),grid[,1])
    grid <- grid[o,]
    evs <- apply(grid,1, function(x){
		 exp(-t[1]*x[1] -t[2]*x[2])
		 })
}

#' Calculates basis vector for 2d cosine transofrm
#' @param n,m number of points in x and y direction
#' @param k vector of length two, giving the index in x and y
#' 	direction
#' @param return.as.vector bool if true, a vector is returned
#' 	if false, a matrix of size nxm is returned
discrete.cosine.2d <- function(n,m, k=c(1,1),
			       return.as.vector=FALSE){
    single.term <- function(x,y,k,n,m){cos((x+0.5)*k[1]*pi/n)*
	    cos((y+0.5)*k[2]*pi/m)}
    mat <- outer(0:(n-1),0:(m-1), single.term, k, n, m)

    if(return.as.vector) mat <- c(mat)
    return(mat)
}

#' centers a matrix
centerscale <- function(x){
    y <- apply(x, 1, function(x)(x - mean(x))/sd(x))
    y <- t(y)
    y[apply(x,1,sd)==0,] <- 0
    return(y)
}
