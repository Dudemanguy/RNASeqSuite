#RNAseq filter function
print("R function to filter outlier data between replicates based on the standard deviation of the difference between normalized vectors")
RNAfilter <- function(x, y){
    if (class(x) != "data.frame") {
        print("Error, please enter a data frame for x")
        }
    if (class(y) != "numeric") {
        print("Error, please enter a numeric value for y")
        }
    else {
        x_filter <- x[rowSums(x) != 0,]
        x_filter[,1] <- x_filter[,1]/norm(data.matrix(x[,1]), type="f")
        x_filter[,2] <- x_filter[,2]/norm(data.matrix(x[,2]), type="f")
        x_diff <- x_filter[,1, drop=FALSE] - x_filter[,2, drop=FALSE]
        mean_x <- mean(x_diff[,1])
        std_x <- sd(x_diff[,1])
        x_compressed <- x_diff[(mean_x - (y*std_x)) <= x_filter[,1, drop=FALSE] & x_filter[,1, drop=FALSE] <= (mean_x + (y*std_x)), 1, drop=FALSE]
        print(paste(nrow(x)-nrow(x_filter), "zero elements discarded"))
        print(paste(nrow(x_diff)-nrow(x_compressed), "outliers removed"))
        return(x_compressed)
        }
}

