plot_dsc <- function(r){
    args <- do.call(rbind, strsplit(r$scenario,"_")) 
    args <- as.data.frame(args)
    names(args) <- c("group", "nDemes", "xscale", "yscale", "s", "angle", "spread")
    r <- cbind(r, args)
    return(r)
}
