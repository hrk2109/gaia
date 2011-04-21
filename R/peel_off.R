`peel_off` <-
function(pvalues_list, threshold, chromosomes, aberrations, discontinuity, hom_threshold){

regions_list <- list();

for (i in 1:length(aberrations)){
	region_chromosome_list <- list();
	cnv_index <- aberrations[i];
	for (j in 1:length(chromosomes) ){		
		property_list <- list();
		chromosome_index <- chromosomes[j];			
		message(".", appendLF = FALSE);	
		curr_pvalue <- pvalues_list[[cnv_index]][[chromosome_index]];
		qvals <- qvalue(curr_pvalue);
		curr_qvalue <- qvals$qvalues;
		tmp_qvalue <- curr_qvalue;			
		tmp_start <- c();
		tmp_end <- c();
		start <- c();
		end <- c();
		qval <- c();
		tmp_list <- list();
		while ( min(tmp_qvalue) < threshold ){

			# Extract the minimum peak of the chromosome
			tmp_list <- search_peaks_in_regions(tmp_qvalue, 1, length(tmp_qvalue), discontinuity[[chromosome_index]], hom_threshold, threshold);
			if (tmp_list[[1]] == -1){
				break;
			}else{
				tmp_start <- tmp_list[[1]];
				tmp_end <- tmp_list[[2]];
				for (k in 1:length(tmp_start)){
					start <- c(start, tmp_start[k]);
					end <- c(end, tmp_end[k]);
					qval <- c(qval, curr_qvalue[tmp_start[k]]);
					curr_pvalue[tmp_start[k]:tmp_end[k]] <- 1;
				}	
			}
			qvals <- qvalue(curr_pvalue);
			tmp_qvalue <- qvals$qvalues;
				
		}
		if (length(start)>0){
			property_list[[1]] <- start;
			property_list[[2]] <- end;
			property_list[[3]] <- qval;
		}else{
			property_list[[1]] <- -1;
			property_list[[2]] <- -1;
			property_list[[3]] <- -1;
		}	
		region_chromosome_list[[chromosome_index]] <- property_list;
	}
	regions_list[[cnv_index]] <- region_chromosome_list;
}
message("\nDone");
return(regions_list);	

}

