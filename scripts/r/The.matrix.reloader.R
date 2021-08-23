the.matrix.reloader <- function(master.data, ID.column.1, ID.column.2=0, species.column, split.by.column=0, split.by.var=0,  OTU.matrix, collapse_top_species = TRUE){
  #Function by Dave Bennett (davidrbennett88@gmail.com)
  #OTU.matrix = your matrix of interactions
  #master.data = your field data, containing all the following information in columns:
  #ID.column.1 = the identifiers that are found BOTH as the first row of OTU.matrix AND as identifiers for each individual in your master.data
  #ID.column.2 = Same as above, this is to allow for when you have multiple samples per individual. Defaults to false
  #species.column = the column containing your species names. Defaults to false for when you don't want to split
  #split.by.column = the column containing the values you want to split OTU.matrix by. In my case, this will be individual study sites. Defaults to false for when you don't want to split
  #split.by.var = the variable found within split.by.column that you want to split by, e.g. a specific site location
  #collapse_top_species = if FALSE, the matrix returned will have a node (column) for each individual bat. If TRUE, matrix will have a node (column) per bat species. Defaults to TRUE
  #print(split.by.var)
  if(split.by.column == 0 & split.by.var !=0){stop("Sorry, split.by.column and split.by.var must either both have values or both be left blank")}
  if(split.by.column != 0 & split.by.var ==0){stop("Sorry, split.by.column and split.by.var must either both have values or both be left blank")}
  if(split.by.column != 0 & split.by.var !=0){  
    if((split.by.var %in% master.data[,split.by.column])==FALSE ){stop("Error: split.by.var not present in split.by.column")}}
  
  #Make lists of the IDs of the sample identifiers and higher-level species
  faeces_species_1 <- as.character(master.data[,ID.column.1])
  names(faeces_species_1) <- master.data[,species.column]
  #cat('faeces_species_1 is ',faeces_species_1, '\n')
  #cat('names(faeces_species_1) is ',names(faeces_species_1), '\n')
  #cat('length(faeces_species_1) is ',length(faeces_species_1), '\n')
  #cat('length(names(faeces_species_1)) is ',length(names(faeces_species_1)), '\n')
  {if(ID.column.2 != 0){
    faeces_species_2 <- as.character(master.data[,ID.column.2])
    names(faeces_species_2) <- master.data[,species.column]}
    else{faeces_species_2<- rep(0, (nrow(master.data)))}
  }
  #return(faeces_species_2)
  
  site_vector <- as.vector(c(0))
  
  if(split.by.var != 0){
    for (myindex in 1:nrow(master.data)) 
    {
      if (master.data[myindex,split.by.column]==split.by.var)
      {site_vector[(myindex*2)] <- as.character(master.data[myindex,ID.column.1])
      
      if(ID.column.2 !=0) {site_vector[(myindex*2)+1] <- as.character(master.data[myindex,ID.column.2])}
      }
      
    }}#This loop makes a vector containing the sample IDs of samples from your desired location, if you have specified that you want the matrix to be split 
  #print(site_vector)
  temp_matrix_1 <- matrix(nrow = nrow(OTU.matrix), ncol=1000)
  #return(temp_matrix_1)
  rownames(temp_matrix_1) <-rownames(OTU.matrix)
  myindex3 <- 0
  {if (split.by.column != 0){
    for (myindex2 in 1:ncol(OTU.matrix)) 
    {
      if ((OTU.matrix[1,myindex2]) %in% site_vector){
        myindex3 <- myindex3+1  
        temp_matrix_1[,myindex3] <- OTU.matrix[,myindex2]
      } #This makes a matrix containing only the samples from the vector made in the loop above (if you specified that you want the matrix to be split)
      myindex2 <- myindex2+1}
    temp_matrix_1<-temp_matrix_1[,1:unique(myindex3)]
  }
    else{temp_matrix_1 <- OTU.matrix[,1:ncol(OTU.matrix)]}
  }
  #print(faeces_species_1)
  #return(temp_matrix_1)
  #print(sum(as.numeric(as.matrix(temp_matrix_1[2:nrow(temp_matrix_1),]))))
  temp_matrix_2 <- temp_matrix_1
  #print(faeces_species_1)
  #for(i in 1:ncol(temp_matrix_2)){
  #  print(temp_matrix_2[1,i])
  #  print(which(faeces_species_1 == temp_matrix_2[1,i]))
  #}
  #print(faeces_species_2)
  #print(temp_matrix_2[1,])
  #Now to make a loop where each of the column heads is given a species name
  for (myindex4 in 1:ncol(temp_matrix_2)){
    if (temp_matrix_2[1,myindex4] %in% faeces_species_1){
      temp_char <- temp_matrix_2[1,myindex4]
      #print(temp_char)
      #print(names((which(faeces_species_1==temp_char))))
      #print('a')
      #print(temp_char)
      sp_name <- names((which(faeces_species_1==temp_char)))
      #print(sp_name)
      temp_matrix_2[1,myindex4] <- sp_name
      #print(myindex4)
      }
    else if (temp_matrix_2[1,myindex4] %in% faeces_species_2){
      #print('b')
      temp_char <- temp_matrix_2[1,myindex4]
      #print(temp_char)
      sp_name <- names((which(faeces_species_2==temp_char)))
      #print(sp_name)
      temp_matrix_2[1,myindex4] <- sp_name
      }
  }
  # make the entries binary rather than weighted
  temp_matrix_2[2:nrow(temp_matrix_2),] <- ifelse(temp_matrix_2[2:nrow(temp_matrix_2),] == 0, 0, 1)
  #return(temp_matrix_2)
  myindex6 <- 1
  if(collapse_top_species==FALSE){
    colnames(temp_matrix_2) <- temp_matrix_2[1,]
    temp_matrix_2 <- temp_matrix_2[-1,]
    temp_matrix_2 <- temp_matrix_2[,-1]
    class(temp_matrix_2) <- 'numeric'
    rownames(temp_matrix_2) <- (OTU.matrix)[-1,1]
    return(temp_matrix_2)
  }
  #print(dim(temp_matrix_2))
  #print(temp_matrix_2[1,])
  temp_matrix_3 <- matrix(nrow = nrow(temp_matrix_2), ncol=ncol(temp_matrix_2))
  rownames(temp_matrix_3) <- rownames(OTU.matrix)
  #print(dim(temp_matrix_3))
  for(myindex5 in 1:ncol(temp_matrix_2))
  {{
    if (temp_matrix_2[1,myindex5] %in% temp_matrix_3[1,]){
      colno <- which(temp_matrix_3[1,]==temp_matrix_2[1,myindex5])
      temp_matrix_3[2:nrow(temp_matrix_3),colno] <- (as.numeric((temp_matrix_3[2:nrow(temp_matrix_3),(which(temp_matrix_3[1,]==temp_matrix_2[1,myindex5])[1])])) 
                                                     + as.numeric((temp_matrix_2[2:nrow(temp_matrix_3),myindex5])))
      
    } else {temp_matrix_3[,myindex6] <- temp_matrix_2[,myindex5] 
    
    myindex6 <- myindex6+1
    }
    
  }
  }#This loop then merges all of the columns by species
  #return(temp_matrix_3)
  #print(temp_matrix_3[1,])
  #Now for the final formatting of the matrix
  temp_matrix_3<- temp_matrix_3[,1:(length(unique(temp_matrix_3[1,])))-1]#its minus 1 as otherwise it'll include an NA in there
  #return(temp_matrix_3)
  
  penultimate_matrix<-data.matrix(temp_matrix_3[2:nrow(temp_matrix_3),])
  colnames(penultimate_matrix)<-temp_matrix_3[1,]
  #print(colnames(penultimate_matrix))
  #return(penultimate_matrix)
  rownames(penultimate_matrix) <-OTU.matrix[2:nrow(OTU.matrix),1]
  #return(penultimate_matrix)
  #print(rownames(penultimate_matrix))
  #return(penultimate_matrix)
  output_matrix <- matrix(as.numeric(unlist(penultimate_matrix[,1:ncol(penultimate_matrix)])),nrow=nrow(penultimate_matrix))
  colnames(output_matrix) <- colnames(penultimate_matrix)#[1:ncol(penultimate_matrix)]
  rownames(output_matrix) <- OTU.matrix[2:nrow(OTU.matrix),1]
  if(0 %in% rowSums(output_matrix)){
    output_matrix <- output_matrix[-which(rowSums(output_matrix)==0),]
  }
  
  return(output_matrix)
}
