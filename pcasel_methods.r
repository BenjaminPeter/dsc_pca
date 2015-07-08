require(lfa)
require(Matrix)
require(MASS)
require(ape)

method.pc1 <- function(input, args){
    pca <- prcomp(t(input))
    input <- t(apply(input, 1, function(x){
                   (x-mean(x)/sd(x))}))

    output <- list()
    output$score <- apply(input, 1, function(x){
                          x %*% pca$x[,1] })
    output$rank <- rank(output$score)
    return(output)
}


method.dufouret.rho <- function(input, args){
    i <- args[['index']] 
    x <- input$svd
    output <- list()
    output$score <- (x$u[,i] * x$d[i])^2
    output$rank <- rank(-output$score)
    return(output)
}
method.dufouret.rho1 <- function(input, args){
    args[['index']] <- 1
    method.dufouret.rho(input, args)
}
method.dufouret.rho2 <- function(input, args){
    args[['index']] <- 2
    method.dufouret.rho(input, args)
}
method.dufouret.h <- function(input, args){
    k <- args$K
    x <- input$svd
    output <- list()
    output$score <- colSums( (t(x$u[,1:k]) * x$d[1:k]) ^ 2)
    output$rank <- rank(-output$score)
    return(output)
}
method.dufouret.hprime <- function(input, args){
    k <- args$K
    x <- input$svd
    output <- list()
    output$score <- rowSums( x$u[,1:k] ^ 2)
    output$rank <- rank(-output$score)
    return(output)
}

method.lfa <- function(input, args){
    k <- args$K
    input <- input$raw
    n <- nrow(input)
    lfa.basis <- lfa(input, k)
    #af.est <- af(input, lfa.basis)

    dev = t(apply(input,1,function(y){
	y1 = as.numeric((y==1) | (y==2))
	y2 = as.numeric(y==2)
	y=c(y1,y2)
	g <- glm(cbind(y, 2-y) ~ -1 + cbind(rep(lfa.basis[,1],2)),
	    family="binomial")
	return(c(g$null.deviance- g$deviance))
      }))

    #ll = function(snp, af){
#	    -sum(snp*log(af) + (2-snp)*log(1-af))
#    }
#    ll.est = sapply(1:nrow(af.est), function(i) {
#		     ll(input[i,], af.est[i,])})
#    return(2 * (null.dev - ll.est))
    output <- list()
    output$score <- dev[1,]
    output$rank <- rank(-output$score)
    return( output)
}

method.gwish <- function(input,args=NULL){
    Sigma.inv <- ginv(cov(input$raw))
    snps <- centerscale(input$raw)
    return(ll.wish(Sigma.inv, snps))
}


method.gwish.tree <- function(input,args=NULL){
    dist <- sim.to.dis(cov(input$raw))
    tree <- nj(dist)
    dist2 <- cophenetic.phylo(tree)
    Sigma.inv <- ginv(dist2)
    snps <- centerscale(input$raw)
    return(ll.wish(Sigma.inv, snps))
}


ll.wish <- function(Sigma.inv,snps){
    A <- a.proj.mat(Sigma.inv)
    rA <- rankMatrix(A)[1]
    WQ <- Sigma.inv %*% A
    rm(A)
    log.Det.WQ <- log.Det(WQ)

    f <- function(snp, log.Det.WQ, WQ, rA){
        ll <- 1/2 * log.Det.WQ +
            rA/4 * log(abs(trace_outer(WQ, snp)))
    }
    ll <-apply(snps, 1, f, log.Det.WQ, WQ, rA)
    ll[is.infinite(ll)] <- NA

    output <- list()
    output$score <- ll
    output$rank <- rank(output$score)
    return(output)
}

method.lfa.gwish <- function(input, args=NULL){
    k <- args$K
    n <- nrow(input$raw)
    lfa.basis <- lfa(input$raw, k)
    snps <- af(input$raw, lfa.basis)
    Sigma.inv <- ginv(cov(input$raw))
    return(ll.wish(Sigma.inv, snps))

}

method.normal <- function(input, args=NULL){
    snps <- centerscale(input$raw)
    Sigma <- cov(snps)
    ll <- dmvnorm(snps, sigma=Sigma, log=T)
    sd0 <- apply(snps,1,sd) == 0
    ll[sd0] <- NA


    output <- list()
    output$score <- ll
    output$rank <- rank(-output$score)
    return(output)
}

method.chisq <- function(input, args=NULL){
    Sigma <- cov(input$raw)
    ll <- dmvnorm(snps, sigma=Sigma, log=T)
    sd0 <- apply(snps,1,sd) == 0
    ll[sd0] <- NA


    output <- list()
    output$score <- ll
    output$rank <- rank(-output$score)
    return(output)
}

method.fst <- function(input, args=NULL){
    mean.freqs <- apply(input$raw, 1, mean)/2
    pi.within <- apply(input$raw, 1, function(snp)mean(snp==1))
    pi.between <- 2 * mean.freqs * (1-mean.freqs)
    output <- list()
    output$score <- 1 - pi.within / pi.between
    output$rank <- rank(-output$score)
    return(output)
}

fst.1locus.2pop <- function(p1, p2, n1, n2){
        alpha1 = 2. * p1 - 2. * p1 * p1
        alpha2 = 2. * p2 - 2. * p2 * p2    
        num = ( p1 -p2)**2 - (n1+n2)*(n1*alpha1 + n2 * alpha2)/(4.*n1*n2*(n1+n2-.1))
        denom = (p1 -p2)**2 + (4.0*n1*n2-n1-n2)*(n1*alpha1 + n2 * alpha2)/(4.*n1*n2*(n1+n2-1.))
        return(num/denom)
}


#' centers a matrix
center <- function(x){
    y <- apply(x, 1, function(x)x - mean(x))
    return(t(y))
}
#' centers a matrix
centerscale <- function(x){
    y <- apply(x, 1, function(x)(x - mean(x))/sd(x))
    y <- t(y)
    y[apply(x,1,sd)==0,] <- 0
    return(y)
}


require(Rcpp)                                                            
                                                                         
cppFunction('double trace_outer(NumericMatrix WQ, NumericVector snp) {
            unsigned int n = snp.size();                              
            double out = 0.;                                          
            double d;                                                 
                                                                      
            for (unsigned int i = 0; i < n; i++){                     
                for (unsigned int j = 0; j < n; j++){                 
                    d = snp(i) - snp(j);                              
                    out += WQ(i,j) * d * d;                           
                }                                                     
            }                                                         
            return out;                                               
}')                                                                   


trace.product<- function(m1, m2) sum(m1 * t(m2))
log.Det <- function(mat, tol=1e-10){
    ev <- Re(eigen(mat)$values)
    lD <- sum(log(ev[ev>tol]))
    return(lD)
}
#' Orthogonal projection matrix from Hank & Hooten 2013
#'
#' Calculates the orthogonal projection for the generalized wishart likelihood
#' @param sigma.inv inverse of the wishart sigma parameter
a.proj.mat <- function(sigma.inv){
    n <- nrow(sigma.inv)
    ones <- cbind(rep(1,n))
    A0 <- (ones %*% t(ones) %*% sigma.inv) / c(t(ones) %*% sigma.inv %*% ones )
    return( diag(c(ones)) - A0 )
}

 ones <- function(n){
        return(rep(1,n))
}


sim.to.dis <- function(q){
        n <- nrow(q)
    -2*q + diag(q) %*% t(ones(n)) + t(t(ones(n))) %*% diag(q)
}

