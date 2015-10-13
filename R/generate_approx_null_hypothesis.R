`generate_approx_null_hypothesis` <-function(obs_data, num_iterations){

# Given the bin_size=1 the number of bin is equal to the number of sample + 1 
num_of_bins <- nrow(obs_data)+1;

# Generate and return the null hypothesis
null_hypothesis <- 0*c(1:num_of_bins);
freq <- 0*c(1:ncol(obs_data));


theta <- apply(obs_data,1,sum);
theta <- theta/ncol(obs_data);

perm_matrix <- matrix(nrow=nrow(obs_data), ncol=num_iterations);

for(i in 1: nrow(perm_matrix)){

	perm_matrix[i,] <- rbinom(ncol(perm_matrix), 1, theta[i])

}

freq <- apply(perm_matrix,2,sum);


uniq_freq <- unique(freq);
for(i in 1:length(uniq_freq)){
	null_hypothesis[uniq_freq[i]+1] <- null_hypothesis[uniq_freq[i]+1]+sum(freq==uniq_freq[i]);
}


null_hypothesis <- null_hypothesis/sum(null_hypothesis);
return(null_hypothesis);

}

