

bad_index <- which(duplicated(cat.matrix(coords_matrix)) == TRUE)

coords_matrix[bad_index, ]

which(duplicated(cat.matrix(as.data.frame(coords_matrix))) == TRUE)

set.seed(155)
temp_coords <- as.data.frame(coords_matrix)
temp_coords$X <- as.vector(coords_matrix[,1]) + rnorm( length(coords_matrix[,1]), 0, 0.000001)
#temp_coords$Y <- as.vector(coords_matrix[,2]) + rnorm(length(coords_matrix[,2]) , 0, 0.0001)

any(duplicated(cat.matrix(temp_coords))) 
which(duplicated(cat.matrix(temp_coords)) == TRUE)



1. figure out what is going wrong with the data that says its duplicated
2. run it the whole code with those problem children
3. try adding tiny amounts of noise to shift the problem data, to see if works
4. Krig.transform.XY look into the fuction to figure it
5. Krig.replicates that 


1. soutir mean meatrix+  + spatial(errors)
2. inverse martix
3. lamba? signal to noise variace, process to nugget variace 
4. nugget variance 
universal krigging
5. 