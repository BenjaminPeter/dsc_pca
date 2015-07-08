#' Add a set of scenarios to a dsc
#' 
#' This function is a wrapper around addScenario, in order
#' to simplify adding a series of simulations with the 
#' same call structure, but different args
#' @param args Arguments that should be passed to
#' 	all simulations
#' @param flexible.args data.frame of argument combinations;
#' 	each combinations corresponds to one scenario
#' @param name function or character or NULL
#'      if function, this will get passed a set of flexible
#' 	args and should return a scenario name
#' 	if NULL, a default title is generated
#' 	if of type character, parameters are added with
#' 	underscore separaters
#' @param ... Arguments passed directly to addScenario
#' @examples 
#' \dontrun{
#' dsc = new.dsc("example_dsc", "example_dsc")
#' datamaker <- function(args){
#' output <- list()
#' output$input <- rnorm(100,args$mean, args$sd)
#' output$meta <- c()
#' }
#' addScenarioGroup(dsc, list(), expand.grid(mean=-2:2, sd=1:3),
#' name="norm", fn=datamaker, seed=1:5)
#' }
#' @export
addScenarioGroup <- function(dsc, args, flexible.args, name=NULL, ...){

    if(is.null(name)) name <- "ScenarioGroup"
    if(is.character(name)){
	name.fn <- function(pars){
	    s <- paste(c(name, pars), collapse='_')
	    return(s)
	}
    }
    if(is.function(name)){
	name.fn <- name
    }

    single.run <- function(params){
	par.names <- names(params)

        for(i in 1:length(params)){
            n <- par.names[i]
            p <- params[i]
            args[[n]] <- p
        }

        addScenario(dsc, name=name.fn(params), args=args, ...)
    }

    apply(flexible.args, 1, single.run)
     
    return()

}

