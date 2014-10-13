`load_cnv` <-
function(segmentation_matrix, markers_list, num_of_samples){

message("Loading Copy Number Data");

# Detect the chromosomes
chromosomes <- as.numeric(sort(unique(names(markers_list))));
chromosomes <- chromosomes[which(!is.na(chromosomes))];

# Detect the aberration kinds
aberration_kinds <- 1:length(unique(segmentation_matrix[,6]));
names(aberration_kinds) <- sort(unique(segmentation_matrix[,6]));

# Detect the samples
samples <- 1:num_of_samples;
if(!is.numeric(segmentation_matrix[,1])){
	sample_names <- unique(segmentation_matrix[,1]);
	for(i in 1:length(sample_names)){
		segmentation_matrix[which(segmentation_matrix[,1]==sample_names[i]),1] <- i;
	}
}

region_list <- list();
# Create the final structure list of the returned list
for(k in 1:length(aberration_kinds)){
	region_list[[ aberration_kinds[k] ]] <- list();
	
	for(i in 1:length(chromosomes) ){
		region_list[[ aberration_kinds[k] ]][[ chromosomes[i] ]] <- matrix(0, length(samples), ncol(markers_list[[ chromosomes[i] ]]));
		rownames(region_list[[ aberration_kinds[k] ]][[ chromosomes[i] ]]) <- samples;
		colnames(region_list[[ aberration_kinds[k] ]][[ chromosomes[i] ]]) <- c(1:ncol(markers_list[[ chromosomes[i] ]]));
	}
}

for(k in 1:length(aberration_kinds)){
	
	ab_ids <- which(segmentation_matrix[,6]==names(aberration_kinds[k]));

	# In this matrix all regions aberrant as aberration_kinds[k] are stored
	tmp_matrix1 <- segmentation_matrix[ab_ids,];
	if(class(tmp_matrix1)=="numeric"){
		tmp_matrix1 <- t(as.matrix(tmp_matrix1));
	}
	for(i in 1:length(chromosomes) ){
		
		# In this matrix all regions aberrant as aberration_kinds[k] for the i-th chromsome are stored
		tmp_matrix2 <- tmp_matrix1[which(tmp_matrix1[,2]==chromosomes[i]),];
		if(class(tmp_matrix2)=="numeric"){
			tmp_matrix2 <- t(as.matrix(tmp_matrix2));
		}
		message(".", appendLF = FALSE);
		for(j in 1:length(samples)){
		 #In this matrix all regions aberrant as aberration_kinds[k] for the i-th chromsome and for the j-th sample are stored
			tmp_matrix3 <- tmp_matrix2[which(tmp_matrix2[,1]==samples[j]),];
			if(class(tmp_matrix3)=="numeric"){
				tmp_matrix3 <- t(as.matrix(tmp_matrix3));
			}
			if(nrow(tmp_matrix3)>0){
				for(t in 1:nrow(tmp_matrix3)){
					start_prob <- tmp_matrix3[t,3];
					end_prob <- tmp_matrix3[t,4];
					
					start_index <- which(markers_list[[ chromosomes[i] ]][1,] == start_prob);
					end_index <- which(markers_list[[ chromosomes[i] ]][2,] == end_prob);
					region_list[[ aberration_kinds[k] ]][[ chromosomes[i] ]][samples[j], start_index:end_index] <- 1;
				}
			}
		}
	}
}	
message("\nDone");
names(region_list) <- names(aberration_kinds);
return(region_list);

}

