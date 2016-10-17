`generate_null_hypothesis` <-
function(obs_data, num_iterations){

# Given the bin_size=1 the number of bin is equal to the number of sample + 1 
num_of_bins <- nrow(obs_data)+1;

# Generate and return the null hypothesis
null_hypothesis <- 0*c(1:num_of_bins);
freq <- 0*c(1:ncol(obs_data));

for (p in 1:num_iterations){
	
	# Apply a Row random pwermutation
	obs_data <- t(apply(obs_data, 1, sample));
	# Compute the frequency obtained by the row random permutation
	freq <- apply(obs_data,2,sum);

	# Create the histogram of the array 'freq' 
	uniq_freq <- unique(freq);
	for(i in 1:length(uniq_freq)){
		null_hypothesis[uniq_freq[i]+1] <- null_hypothesis[uniq_freq[i]+1]+sum(freq==uniq_freq[i]);
	}
}

null_hypothesis <- null_hypothesis/sum(null_hypothesis);
return(null_hypothesis);

}

