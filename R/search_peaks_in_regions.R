`search_peaks_in_regions` <-
function(qvalues, start_region, end_region){

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
list_peaks[[1]] <- start;
list_peaks[[2]] <- end;
return(list_peaks);	

}

