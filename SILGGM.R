####Harshini####
###This Program is for scaled LASSO
###LASSO regression. The file gets saved in Cytoscape format
### File can be imported into Cytoscape and generate network

setwd('C:/Users/harsh/Desktop/205')

###SILGGM: Package for LASSO regression
library(SILGGM)

##Read gene expression data
MyData <- read.csv(file="WS_Light_PCA.csv",header=TRUE)

##Transpose data. Columns are genes. Rows are expression values.
transpose_data <- t(MyData)
## Column names (headers of the file) must be genes.
colnames(transpose_data)=transpose_data[1,]
transpose_data=transpose_data[-1,]

##Save the transposed data into a table. Cannot use .csv format as the column limit is only 16324
##If you have large data like 20,000 samples, it is better to use only table format
write.table(transpose_data, file = "data")

### Read the saved table
df<-read.table("data")

##Convert that table into a matrix. SILGGM uses only matrix format
mat <- as.matrix(df)

##Here the method used is Scaled LASSO. There are 5 methods in total
## You can refer the documentation.
##The output file generated will be saved in the Cytoscape formatJ
Pass_grn <- SILGGM(mat,method = "GFC_SL", alpha = 1, cytoscape_format = TRUE, csv_save = TRUE)


