filter_mutations <- function(a_matrix, depth, vaf, strand_bias = T){
	a_filter <- a_matrix[,"total"] > depth &
		a_matrix[,"vaf"] > vaf
	if(strand_bias){
		a_filter <- a_filter & !(a_matrix[,"strand_bias"] < 0.05)
	}
	return(a_matrix[a_filter,])
}

filter_mutations_for_eqtl <- function(a_matrix, depth, strand_bias = T){
	a_filter <- a_matrix[,"total"] > depth 
	if(strand_bias){
		a_filter <- a_filter & !(a_matrix[,"strand_bias"] < 0.05)
	}
	return(a_matrix[a_filter,])
}