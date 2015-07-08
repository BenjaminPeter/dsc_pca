
score.top100 <- function(...){
    return(score.topn(..., n=100))
}
score.top500 <- function(...){
    return(score.topn(..., n=500))
}
score.top50 <- function(...){
    return(score.topn(..., n=50))
}
score.topn <- function(data, output, n=100, reverse=F){
    if(reverse){
	r <- length(output$score) - output$score
	positives <- r < n
    } else{
	positives <- output$rank < n
    }
    true.posiitves <- sum(data$meta$issel & positives)
    false.posiitves <- sum((!data$meta$issel) & positives)
    return(list(true.posiitves=true.posiitves,
                false.posiitves=false.posiitves))
}

score.bot100 <- function(...){
    return(score.topn(..., n=100, reverse=T))
}
