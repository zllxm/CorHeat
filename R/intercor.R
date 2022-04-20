
#' Computes a matrix of Spearman's rho or Pearson's r rank correlation coefficients for all possible pairs of columns of a matrix, filter based on the threshold of P and r cutoff, and output the matrix.
#' @param x a numeric matrix
#' @param y a numeric matrix
#' @param method specifies the type of correlations to compute
#' @param pcutoff Specifies a threshold for the p value
#' @param rcutoff Specifies a threshold for the r value
#' @importFrom psych corr.test
#' @importFrom dplyr intersect
#' @export
intercor<-function(x,y,method,pcutoff,rcutoff){
    cor_result<-psych::corr.test(x,y,method= method)
    r_result<-cor_result$r
    p_result<-cor_result$p
    p_col_row <-data.frame(which(p_result[] < pcutoff, arr.ind = TRUE))
    r_col_row <-data.frame(which(abs(r_result[]) > rcutoff, arr.ind = TRUE))
    r_p_inter<-dplyr::intersect(p_col_row,r_col_row)
    row_n<-unique(r_p_inter$row)
    col_n<-unique(r_p_inter$col)
    r_filter_result<- r_result[row_n,col_n]
    p_filter_result<- p_result[row_n,col_n]
    filter_result<- list(r_filter_result=r_filter_result,p_filter_result=p_filter_result)
    return(filter_result)
  }

