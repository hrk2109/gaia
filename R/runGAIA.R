`runGAIA` <-
function(cnv_obj, markers_obj, output_file_name="",  aberrations=-1, chromosomes=-1, num_iterations=10, threshold=0.25, hom_threshold=0.12, approximation=FALSE){

message("\nPerforming Data Preprocessing");
# Chromosome indexing array
if(chromosomes==-1){ 
	# We apply the algorithm to each chromosome
	chromosomes <- as.numeric(sort(unique(names(markers_obj))));
	chromosomes <- chromosomes[which(!is.na(chromosomes))];

}else{
	known_chr <- as.numeric(names(markers_obj));
	if( (length(known_chr[chromosomes])!= length(chromosomes)) || (sum(is.na(known_chr[chromosomes]))>0) ){
		error_string <- "Error in the list of chromosomes passed as argument.\n";
		error_string <- cat(error_string, "The list of the chromosomes that can be analyzed follows:\n", known_chr, "\n");
		stop(error_string, call.=FALSE);
	}
}


# Aberration Kind indexing array
if(aberrations==-1){ 
	# We apply the algorithm to each aberrations
	aberrations <- as.numeric(names(cnv_obj));
	names(aberrations) <- names(cnv_obj);
	if(length(aberrations)>0 && aberrations[1]==0){
		aberrations <- aberrations+1;
	}
	
}else{
	aberrations <- as.numeric(names(cnv_obj));
	names(aberrations) <- names(cnv_obj);
	if(length(aberrations)>0 && aberrations[1]==0){
		aberrations <- aberrations+1;
	}
	known_aberr <- as.numeric(names(cnv_obj));
	if( (length(known_aberr[aberrations])!= length(aberrations)) || (sum(is.na(known_aberr[aberrations])>0)) ){
		error_string <- "Error in the list of aberrations passed as argument.\n";
		error_string <- cat(error_string, "The aberrations that can be analyzed follow:\n", known_aberr, "\n");
		stop(error_string, call.=FALSE);
	}
}

# Discontinuity Matrix for Homogeneous peel-off
discontinuity <- list();
if(length(aberrations)==2 && hom_threshold>=0){
	message("\nDone");
	message("\nComputing Discontinuity Matrix");
	for (i in 1:length(chromosomes)){
		message(".", appendLF = FALSE);
		tmp <- cnv_obj[[2]][[chromosomes[i]]]-cnv_obj[[1]][[chromosomes[i]]];
		tmp_vec <- 0*c(1:(ncol(tmp)-1));
	
		for(k in 1:(ncol(tmp)-1)){
			for(z in 1:nrow(tmp)){
				tmp_vec[k] <- tmp_vec[k]+ abs(tmp[z,k]-tmp[z,k+1]);
			}
		}
		discontinuity[[chromosomes[i]]] <- tmp_vec/nrow(cnv_obj[[2]][[chromosomes[i]]]);
	}
	
}else{
	if(length(aberrations)!=2 && hom_threshold>=0){
		message("\nHomogeneous cannot be applied on the data (data must contain exactly two different kinds of aberrations)\n");
	}
	for (i in 1:length(chromosomes)){
		tmp <- cnv_obj[[1]][[chromosomes[i]]]-cnv_obj[[1]][[chromosomes[i]]];
		tmp_vec <- 0*c(1:(ncol(tmp)-1));		
		discontinuity[[chromosomes[i]]] <- tmp_vec;
		hom_threshold = -1;
	}
}
message("\nDone");

# The bin size is fixed to 1
# Generate The Null hypothesis for each chromosome and aberration
null_hypothesis_list <- list();
null_hypothesis_chromosome_list <- list();

message("Computing Probability Distribution");
for (k in 1:length(aberrations)){
	null_hypothesis_chromosome_list <- list();
	for (i in 1:length(chromosomes)){
		
		aberrations_index <- (aberrations[k]);
		chromosome_index <- chromosomes[i];
				
		obs_data <- cnv_obj[[aberrations_index]][[chromosome_index]];
		message(".", appendLF = FALSE);
		if(approximation==FALSE){
			null_hypothesis_chromosome_list[[chromosome_index]] <- generate_null_hypothesis(obs_data, num_iterations);
		}
		if(approximation==TRUE){
			null_hypothesis_chromosome_list[[chromosome_index]] <- generate_approx_null_hypothesis(obs_data, num_iterations);
		}
		
	}	
	null_hypothesis_list[[aberrations_index]] <- null_hypothesis_chromosome_list;
}

# Compute the p-value curve for each chromosome and aberration
pvalue_distribution_list <- list();
for (k in 1:length(aberrations)){
	tmp_pvalue_list <- list();
	for (i in 1:length(chromosomes)){
	message(".", appendLF = FALSE);
	tmp_pvalue_list[[chromosomes[i]]] <- rev(cumsum(rev(null_hypothesis_list[[aberrations[k]]][[chromosomes[i]]])));			
	}
	pvalue_distribution_list[[aberrations[k]]] <- tmp_pvalue_list;
}
message("\nDone");


# Compute the pvalue for each cnv, chromosome and marker
pvalues_list <- list();

message("Assessing the Significance of Observed Data");
for (k in 1:length(aberrations)){
	pvalue_chromosome_list <- list();
	for (i in 1:length(chromosomes)){
		message(".", appendLF = FALSE);
		aberrations_index <- (aberrations[k]);
		
		chromosome_index <- chromosomes[i];
		# the p-value curve for the k-th aberration and i-th chromosome
		curr_pvalue <- pvalue_distribution_list[[aberrations_index]][[chromosome_index]];	

		# Observed data for k-th aberration and i-th chromosome
		curr_obs_data <- cnv_obj[[aberrations_index]][[chromosome_index]];
		
		# Compute the absolute frequency for the current observed data
		obs_freq <- apply(curr_obs_data, 2, sum);

		# Assign to each observed point the respective computed p-value
		pvalue_chromosome_list[[chromosome_index]] <- round(curr_pvalue[obs_freq[]+1],7);
		
	}	
	pvalues_list[[aberrations_index]] <- pvalue_chromosome_list;
}
message("\nDone");
# Generate the .gistic file
if(length(aberrations)==2){
	gistic_label <- c("Del", "Amp");
	message("Writing ", paste(output_file_name, ".igv.gistic", sep=""), " File for Integrative Genomics Viewer (IGV) Tool");
	gistic <- matrix(nrow=0, ncol=8);
	colnames(gistic) <- c("Type", "Chromosome", "Start", "End", "q-value", "G-score", "average amplitude", "frequency")
	for (k in 1:length(aberrations)){
		tmp_pvalue_list <- list();
		for (i in 1:length(chromosomes)){
			message(".", appendLF = FALSE);
			curr_qval <- as.numeric(pvalues_list[[aberrations[k]]][[chromosomes[i]]]);
			curr_qval <- qvalue(curr_qval);
			start <- 1;
			for(z in 2:(length(curr_qval))){
				if(curr_qval[z-1]!=curr_qval[z]){
					end <- z-1;
					print_qval <- curr_qval[z-1];
					
					gistic <- rbind(gistic, c(gistic_label[k], chromosomes[i], markers_obj[[chromosomes[i]]][1,start],  
								markers_obj[[chromosomes[i]]][2,end],print_qval,0,0,0)); 
					start <- z;						
				}
			}
			end <- length(curr_qval);
			print_qval <- curr_qval[end];
			gistic <- rbind(gistic, c(gistic_label[k], chromosomes[i], markers_obj[[chromosomes[i]]][1,start],  								markers_obj[[chromosomes[i]]][2,end], print_qval,0,0,0));
				
			
		}
	}
	write.table(gistic, paste(output_file_name, ".igv.gistic", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, eol="\n", quote=FALSE)
	message("\nDone");
}

# Running the peel-off algorithm
if(hom_threshold>=0){
	message("Running Homogeneous peel-off Algorithm With Significance Threshold of ", threshold,  " and Homogenous Threshold of ", hom_threshold);
}else{
	message("Running Standard peel-off Algorithm With Significance Threshold of ", threshold);
}
significant_regions_list <- peel_off(pvalues_list, threshold, chromosomes, aberrations, discontinuity, hom_threshold);
#return(list(significant_regions_list, chromosomes));
if(output_file_name!=""){
	# Generation of output file
	message("\nWriting Output File \'", output_file_name, "\' Containing the Significant Regions\n", sep="");
	results <- write_significant_regions(markers_obj, significant_regions_list, output_file_name, chromosomes, aberrations);
	message("File \'", output_file_name, "\' Saved\n", sep="");
	
}else{
	results <- write_significant_regions(markers_obj, significant_regions_list, output_file_name, chromosomes, aberrations);
}
return(results);
}

