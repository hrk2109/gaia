`search_peaks_in_regions` <-
function(qvalues, start_region, end_region, discontinuity, hom_threshold, threshold){

index <- 0;
in_min <- FALSE;
start <- c();
end <- c();

list_peaks <- list();

min_region <- min(qvalues[start_region:end_region]);

for ( i in start_region:end_region ){

	if ( in_min == FALSE){
		if ( qvalues[i]==min_region ){
			index <- index +1;
			start[index] <- i;
			in_min <- TRUE;
		}
	}else{
		if ( (qvalues[i] != min_region) && (qvalues[i-1] == min_region) ){
			end[index] <- i-1;
			in_min <- FALSE;
			
		}
	}	
}

if (in_min==TRUE){
	end[index] <- end_region;
}

# Extend boundaries of selected peak
if(hom_threshold>=0){
	for(z in 1:length(end)){
		curr_start <- start[z];
		curr_end <- end[z];

		max_left <-  max(discontinuity[curr_start:start_region]);
		new_start <- curr_start;

		if( max_left >= hom_threshold && min(qvalues[curr_start:start_region])<threshold){
			for(i in curr_start:start_region){
				if(qvalues[i]>threshold){
						break;
				}
				if(discontinuity[i] >= hom_threshold && qvalues[i]<threshold){
					new_start <- i+1;
					break;
				}
			}
		}	
	
		new_end <- curr_end;
		if(curr_end<end_region){
			max_right <-  max(discontinuity[curr_end:(end_region-1)]);

			if(max_right >= hom_threshold && min(qvalues[curr_end:end_region])<threshold){
				for(i in curr_end:(end_region-1)){
					if(qvalues[i]>threshold){
						break;
					}
					if(discontinuity[i] >= hom_threshold && qvalues[i]<threshold){
						new_end <- i;
						break;
					}
				}
			}
		}
	}	
}
list_peaks[[1]] <- start;
list_peaks[[2]] <- end;
return(list_peaks);	
}

