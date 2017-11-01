get_forced_beta <- function(cem, network_type){
	num_samples <- ncol(expr_data(cem))
	num_tests <- c(20, 30, 40, 60)
	if(network_type=="unsigned"){
		beta_values <- c(9, 8, 7, 6)
		if(num_samples < 20){
			beta <- 10
		}else{
			beta <- beta_values[which(num_tests == max(num_tests[num_tests <= num_samples]))]
		}
	}else if(network_type=="signed"){
		beta_values <- c(18, 16, 14, 12)
		if(num_samples < 20){
			beta <- 20
		}else{
			beta <- beta_values[which(num_tests == max(num_tests[num_tests <= num_samples]))]
		}
	}
	return(beta)
}

get_args <- function(cem, vars){
	system_calls <- lapply(sys.calls(), as.character)
	go_back <- ifelse(".local" %in% unlist(system_calls), 2, 1)

	func <- as.character(sys.calls()[[sys.nframe()-go_back]][[1]])
	if(class(get(func)) == "function"){
		command <- deparse(sys.call(-1))
	}else if(class(get(func)) == "nonstandardGenericFunction"){
		command <- deparse(sys.call(-2))
	}
	exceptions <- c("expr", "sample_annot", "gmt", "interactions")
	vars[["cem"]] <- NULL
	vars[["results"]] <- NULL
	vars[names(vars) %in% exceptions] <- "input"
		    
	inputs <- cem@input_params

	inputs[[func]] <- vars
			    
	cem@input_params <- inputs
				    
	calls <- cem@calls
	calls[[length(calls) + 1]] <- command
	if(length(calls) == 1){
		names(calls) <- func
	}else{
		names(calls)[length(calls)] <- func
	}
						    
	cem@calls <- calls
	return(cem)
}
