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
