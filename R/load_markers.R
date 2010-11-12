`load_markers` <-
function(marker_matrix){

message("Loading Marker Informations");

chromosomes <- sort(unique(marker_matrix[,2]));


chromosome_marker_list <- list();

# Check if in the file there is the fourth column for the end position of the probes; 
end_position <- FALSE;
if(ncol(marker_matrix)==4){
	end_position <- TRUE;
}

for( i in 1:length(chromosomes) ){
	chr_ids <- which(marker_matrix[,2]==chromosomes[i]);
	tmp_matrix <- matrix(0,2,length(chr_ids));
	tmp_matrix[1,] <- marker_matrix[chr_ids,3];
	if(end_position){
		tmp_matrix[2,] <- marker_matrix[chr_ids,4];
	}else{
		tmp_matrix[2,] <- tmp_matrix[1,];
	}
	chromosome_marker_list[[chromosomes[i]]] <- tmp_matrix;
	message(".", appendLF = FALSE);
}
message("\nDone");
names(chromosome_marker_list) <- chromosomes;
return(chromosome_marker_list);
}

