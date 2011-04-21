`write_significant_regions` <-
function(markers_list, regions_list, output_file_name, chromosomes, aberrations){

# This function writes into a tab delimited file the significant aberrant regions.

num_of_cnv <- length(aberrations);
num_of_chromosome <- length(chromosomes);

aberration_list <- list();
index <- 0;
aberration_matrix <- matrix(nrow=0, ncol=6);

for (i in 1:num_of_cnv){
	cnv_index <- aberrations[i];
	for (j in 1:num_of_chromosome){
		chromosome_index <- chromosomes[j];

		markers <- markers_list[[chromosome_index]];
		
		sign_regions <- regions_list[[cnv_index]][[chromosome_index]];
		start <- sign_regions[[1]];
		end <- sign_regions[[2]];
		qval <- sign_regions[[3]];
		if (start[1]!=-1){
			for ( k in 1:length(start) ){
				index <- index + 1;		
	region_line <- c(chromosome_index, names(aberrations[i]), markers[1,start[k]], markers[2,end[k]], ((markers[2,end[k]]-markers[1,start[k]])+1), qval[k]);
			aberration_matrix <- rbind(aberration_matrix, region_line);
			}
		}
	}
}
if(nrow(aberration_matrix)>0){
	rownames(aberration_matrix) <- c(1:nrow(aberration_matrix));
}
colnames(aberration_matrix) <- c("Chromosome", "Aberration Kind", "Region Start [bp]", "Region End [bp]", "Region Size [bp]", "q-value");
if(output_file_name!=""){
	write.table(aberration_matrix, output_file_name, quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE);
}
return(aberration_matrix);
}

