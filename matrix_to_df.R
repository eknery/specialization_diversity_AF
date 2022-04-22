### setting function to convert matrix to dataframe

matrix_to_df = function(matrix){
  # AF matrix
  af_matrix = matrix[rownames(matrix) %in% AF, colnames(matrix) %in% AF]
  af_values = c()
  for (i in 1:nrow(af_matrix)){
    for(j in i:ncol(af_matrix)){
      if(i == j){next}
      af_values= c(af_values, af_matrix[i,j])
    }
  }
  
  # WS matrix
  ws_matrix = matrix[!rownames(matrix) %in% AF, !colnames(matrix) %in% AF]
  ws_values = c()
  for (i in 1:nrow(ws_matrix)){
    for(j in i:ncol(ws_matrix)){
      if(i == j){next}
      ws_values= c(ws_values, ws_matrix[i,j])
    }
  }
  
  # organizing into dataframe
  af_name = rep("AF-endemic", length(af_values))
  ws_name = rep("widespread", length(ws_values))
  distribution = c(af_name,ws_name)
  all_values = c(af_values,ws_values)
  df = data.frame(distribution, all_values)
  return(df)
}