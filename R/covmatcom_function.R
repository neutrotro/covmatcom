# ___________________________________________________________________________
#|                                                                           |
#|                                   covmatcom                               | 
#| __________________________________________________________________________| 
#  
#' R function Comparing Two Datasets' Covariance Matrices 
#'

#' This function takes as input two datasets with measurements (rows) for a
#' set of variables (columns; the same number of variables in the two 
#' datasets), calculates the two datasets' covariance matrices, measures the
#' overall differentiation between these matrices and separates the  
#' contribution of matrix orientation and matrix shape to the differentiation.
#' The function produces as output: i) s1, an overall measure of 
#' differentiation, ii) s2, measuring the contribution of differences in 
#' orientation to s1, and iii) s3, measuring the contribution of differences
#' in matrix shape to s1.
                                                                          
#' See:  Garcia C. 2012. A simple procedure for the comparison of covariance 
#'   matrices. BMC Evolutionary Biology 12:1-17 for details. 
#'
#' @param data.a Dataset to calculate the first cov matrix to be compared
#' @param data.b Dataset to calculate the second matrix to be compared
#' @param nci.boo number of bootstrap resamples to obtain upper confidence
#'        intervals. If = 0, no C.I. are calculated. Defaults to  = 0.
#' @keywords covMatrix, comparison, orientation, shape
#' @author Carlos Garcia, \email{carlos.garcia.suarez@usc.es}
#' @references Garcia C. 2012. A simple procedure for the comparison of covariance 
#'   matrices. BMC Evolutionary Biology 12:1-17
#'
#' @export
#' @return A list including "NumberOfVariables", "s_statistics" and
#'                           "Null_Dist_Upper_Quantiles"

#'    NOTE: The function obtains these null distribution quantiles by 
#'  calculating the s statistics for pairs of bootstrap resamples taken from a
#'   single data set. A set of s statistics' quantiles is obtained for dataA 
#'  and another for dataB.



# covmatcom_function()
# _____________________________________________________________________________ 
#|                                                                             | 
#| Carlos Garcia                                                               | 
#| CIBUS Campus Sur                                                            | 
#| Universidade de Santiago de Compostela                                      |
#| 15782 A Corunha                                                             | 
#| Galiza - Spain                                                              | 
#|                                                                             | 
#| carlos.garcia.suarez@usc.es                                                 | 
#|_____________________________________________________________________________| 
                                                                             
# ****************************************************************************** 
# USAGE 

# INPUT:  Two numeric arrays of two dimensions, observations in rows, variables 
  # in columns 

# OUTPUT: 1) "result.out": object including: 
  # - "NumVar": Number of variables 
  # - "Statistics": S1, S2 and S3 statistics values 
  # - "Quantiles": 0.999, 0.990 and 0.950 quantiles for the three statistics 
    # from null distributions obtained by bootsrapping within each data set.
 


# CALL: "result.out<-covmatcom(data.a,data.b,nci.boo=0)" 

  


# ***************************************************************************** 
# ***************************************************************************** 
# Function  
# ............................................................................. 
# ............................................................................. 

covmatcom <- function(data.a,data.b,nci.boo=0){ 
  
  # data.a: first data file to calculate the covariance matrix with, 
     # observations in rows, variables in columns 
  # data.b: second data file to calculate the covariance matrix with, 
     # observations in rows, variables in columns 
  # nci.boo: number of resamples to calculate null distributions'
     # bootstrap confidence intervals for the differentiation statistics
     # if zero, no intervals are calculated. 


# ***************************************************************************** 
# Number of variables checking
ncol.a<-ncol(data.a)
ncol.b<-ncol(data.b)
if(ncol.a ==! ncol.b) {print.default("ERROR: incompatible variable numbers")}


# ***************************************************************************** 
# Statistics' observed values
# ***************************************************************************** 

# .............................................................................
# Storage vectors
# .............................................................................
  dif.a<-numeric(ncol.a)
  dif.b<-numeric(ncol.a) 
  hori<-numeric(ncol.a) 
  diag<-numeric(ncol.a) 
  vexp.aa<-numeric(ncol.a)
  vexp.ab<-numeric(ncol.a)
  vexp.ba<-numeric(ncol.a)
  vexp.bb<-numeric(ncol.a) 

# .............................................................................
# Eigenvector analyses
# .............................................................................
  

nrow.a<-nrow(data.a)
nrow.b<-nrow(data.b)

# Covariance matrices
 mat.a<-cov(data.a) 
 mat.b<-cov(data.b) 

eigen.a<-eigen(mat.a) 
eigen.b<-eigen(mat.b) 

# Eigenvectors
egvc.a<-eigen.a[["vectors"]] 
egvc.b<-eigen.b[["vectors"]] 

# Eigenvalues
egvl.a<-eigen.a[["values"]] 
egvl.b<-eigen.b[["values"]] 
vartot.a<-sum(egvl.a)
vartot.b<-sum(egvl.b) 

# .............................................................................
# Statistics' values 
# .............................................................................

 
  # Individual observations' values for own/reciprocal eigenvectors 
  egvc.aa<-(data.a)%*%egvc.a 
  egvc.ab<-(data.a)%*%egvc.b 
  egvc.ba<-(data.b)%*%egvc.a 
  egvc.bb<-(data.b)%*%egvc.b 

 # Amount of variance associated with every own/reciprocal eigenvector
  # in data samples a and b
  for(ip in 1:ncol.a) { 
  vexp.aa[ip]<-var(egvc.aa[,ip])/vartot.a 
  vexp.ab[ip]<-var(egvc.ab[,ip])/vartot.a 
  vexp.ba[ip]<-var(egvc.ba[,ip])/vartot.b 
  vexp.bb[ip]<-var(egvc.bb[,ip])/vartot.b 
               
# own/reciprocal comparisons of variance explained   
  
  dif.a[ip] <- vexp.ba[ip]-vexp.aa[ip] 
  dif.b[ip] <- vexp.ab[ip]-vexp.bb[ip] 
  hori[ip] <- (vexp.aa[ip]+vexp.ab[ip]) - (vexp.ba[ip]+vexp.bb[ip]) 
  diag[ip] <- ( vexp.aa[ip]+vexp.bb[ip]) - (vexp.ba[ip]+vexp.ab[ip]) 


                       } # end of "for(ip ..." 

s1o<-2*(sum(dif.a^2)+sum(dif.b^2)) 
s2o<-sum(diag^2) 
s3o<-sum(hori^2) 


# All calculations are divided by 8 to obtain 0 to 1 values for s1. 
# This is because the maximum possible values for sum(dif.a²) and sum(dif.b²)
# are = 2. These values would occur if all variance in each sample were
# explained by single, different eigenvectors, so that some, say dif.a[ip]
# would be equal to +1 and another, dif.a[ip'], to -1. Thus the maximum  
# possible value for sum(dif.a²) is = 2 (and also =2 for sum(dif.b²)), 
# and is = 2*(2+2) = 8 for s1o

s1<-s1o/8 # total differentiation statistic
s2<-s2o/8 # orientation contribution statistic
s3<-s3o/8 # shape contribution statistic

s.obs<-c(s1,s2,s3)
names(s.obs)<-c("s1-total","s2-orientation","s3-shape")



  quantiles=NULL
if(nci.boo > 0) {
# ***************************************************************************** 
# Bootstrap Confidence intervals
# ***************************************************************************** 

# .............................................................................
# Storage vectors
# .............................................................................
 dif.a1r<-numeric(ncol.a)
 dif.a2r<-numeric(ncol.a)
 dif.b1r<-numeric(ncol.a)
 dif.b2r<-numeric(ncol.a)
   hori.ar<-numeric(ncol.a)
   hori.br<-numeric(ncol.a)
 diag.ar<-numeric(ncol.a)
 diag.br<-numeric(ncol.a)
   vexp.a1a1r<-numeric(ncol.a)
   vexp.a1a2r<-numeric(ncol.a)
   vexp.a2a1r<-numeric(ncol.a)
   vexp.a2a2r<-numeric(ncol.a) 
   vexp.b1b1r<-numeric(ncol.a)
   vexp.b1b2r<-numeric(ncol.a)
   vexp.b2b1r<-numeric(ncol.a)
   vexp.b2b2r<-numeric(ncol.a) 

s1.ar<-numeric(nci.boo ) 
s2.ar<-numeric(nci.boo )  
s3.ar<-numeric(nci.boo ) 
s1.br<-numeric(nci.boo )  
s2.br<-numeric(nci.boo )  
s3.br<-numeric(nci.boo ) 

# .............................................................................
# Bootstrap resampling loop
# .............................................................................
ref.a<-c(1:nrow.a)
ref.b<-c(1:nrow.b)

print.default("Bootstrapping ...")
 
for(iboo in 1:nci.boo)       { 

ref.a1r<-sample(ref.a,nrow.a,replace=TRUE)
ref.a2r<-sample(ref.a,nrow.a,replace=TRUE)
ref.b1r<-sample(ref.b,nrow.a,replace=TRUE)
ref.b2r<-sample(ref.b,nrow.a,replace=TRUE)

data.a1r<-data.a[ref.a1r,]
data.a2r<-data.a[ref.a2r,]
data.b1r<-data.a[ref.b1r,]
data.b2r<-data.a[ref.b2r,]


# Covariance matrices
 mat.a1r<-cov(data.a1r) 
 mat.a2r<-cov(data.a2r) 
 mat.b1r<-cov(data.b1r) 
 mat.b2r<-cov(data.b2r)

# Eigen analysis
eigen.a1r<-eigen(mat.a1r) 
eigen.a2r<-eigen(mat.a2r) 
eigen.b1r<-eigen(mat.b1r) 
eigen.b2r<-eigen(mat.b2r) 

# Eigenvectors
egvc.a1r<-eigen.a1r[["vectors"]] 
egvc.a2r<-eigen.a2r[["vectors"]] 
egvc.b1r<-eigen.b1r[["vectors"]] 
egvc.b2r<-eigen.b2r[["vectors"]] 

# Eigenvalues
egvl.a1r<-eigen.a1r[["values"]] 
egvl.a2r<-eigen.a2r[["values"]] 
egvl.b1r<-eigen.b1r[["values"]] 
egvl.b2r<-eigen.b2r[["values"]] 
vartot.a1r<-sum(egvl.a1r)
vartot.a2r<-sum(egvl.a2r)
vartot.b1r<-sum(egvl.b1r) 
vartot.b2r<-sum(egvl.b2r)
 

  # Individual observations' values for own/reciprocal eigenvectors 
  egvc.a1a1r<-(data.a1r)%*%egvc.a1r 
  egvc.a1a2r<-(data.a1r)%*%egvc.a2r 
  egvc.a2a1r<-(data.a2r)%*%egvc.a1r 
  egvc.a2a2r<-(data.a2r)%*%egvc.a2r 
  egvc.b1b1r<-(data.b1r)%*%egvc.b1r 
  egvc.b1b2r<-(data.b1r)%*%egvc.b2r 
  egvc.b2b1r<-(data.b2r)%*%egvc.b1r 
  egvc.b2b2r<-(data.b2r)%*%egvc.b2r 

 # Amount of variance associated with every own/reciprocal eigenvector
  # in data samples a and b
  for(ip in 1:ncol.a) { 
  vexp.a1a1r[ip]<-var(egvc.a1a1r[,ip])/vartot.a1r 
  vexp.a1a2r[ip]<-var(egvc.a1a2r[,ip])/vartot.a1r 
  vexp.a2a1r[ip]<-var(egvc.a2a1r[,ip])/vartot.a2r 
  vexp.a2a2r[ip]<-var(egvc.a2a2r[,ip])/vartot.a2r 

  vexp.b1b1r[ip]<-var(egvc.b1b1r[,ip])/vartot.b1r 
  vexp.b1b2r[ip]<-var(egvc.b1b2r[,ip])/vartot.b1r 
  vexp.b2b1r[ip]<-var(egvc.b2b1r[,ip])/vartot.b2r 
  vexp.b2b2r[ip]<-var(egvc.b2b2r[,ip])/vartot.b2r 
 
  dif.a1r[ip] <- vexp.a2a1r[ip]-vexp.a1a1r[ip] 
  dif.a2r[ip] <- vexp.a1a2r[ip]-vexp.a2a2r[ip] 
hori.ar[ip] <- (vexp.a1a1r[ip]+vexp.a1a2r[ip]) -
 (vexp.a2a1r[ip]+vexp.a2a2r[ip]) 
diag.ar[ip] <-( vexp.a1a1r[ip]+vexp.a2a2r[ip]) -
 (vexp.a2a1r[ip]+vexp.a1a2r[ip]) 
 
  dif.b1r[ip] <- vexp.b2b1r[ip]-vexp.b1b1r[ip] 
  dif.b2r[ip] <- vexp.b1b2r[ip]-vexp.b2b2r[ip] 
hori.br[ip] <- (vexp.b1b1r[ip]+vexp.b1b2r[ip]) -
 (vexp.b2b1r[ip]+vexp.b2b2r[ip]) 
diag.br[ip] <-( vexp.b1b1r[ip]+vexp.b2b2r[ip]) -
 (vexp.b2b1r[ip]+vexp.b1b2r[ip]) 
                       } 


s1.apr<-2*(sum(dif.a1r^2)+sum(dif.a2r^2)) 
s2.apr<-sum(diag.ar^2) 
s3.apr<-sum(hori.ar^2) 

s1.bpr<-2*(sum(dif.b1r^2)+sum(dif.b2r^2)) 
s2.bpr<-sum(diag.br^2) 
s3.bpr<-sum(hori.br^2) 

# All calculations are divided by 8 to obtain 0 to 1 values for s1.apr
# and s1.bpr. This is because the maximum possible values for sum(dif.a1r²),
# sum(dif.a2r²), sum(dif.b1r²) and sum(dif.b2r²) are = 2. These values would
# occur if all variance in each sample were explained by single different
# eigenvectors, so that some, say dif.a1r[ip] would be equal to +1 and
# another, dif.a1r[ip'], to -1. Thus the maximum possible value for
# sum(dif.a1r²) is = 2 (and also =2 for sum(dif.a21r²)), 
# and is = 2*(2+2) = 8 for s1apr and s1bpr.


s1.ar[iboo]<-s1.apr/8 # total differentiation statistic
s2.ar[iboo]<-s2.apr/8 # orientation contribution statistic
s3.ar[iboo]<-s3.apr/8 # shape contribution statistic

s1.br[iboo]<-s1.bpr/8 # total differentiation statistic
s2.br[iboo]<-s2.bpr/8 # orientation contribution statistic
s3.br[iboo]<-s3.bpr/8 # shape contribution statistic

  
            }   # "for(iboo..."  end

print.default(" … done.")


# .............................................................................
# Null distributions' quantiles
# .............................................................................

qval<-c(0.999,0.990,0.950)

quant.1a<-quantile(s1.ar,qval)
quant.2a<-quantile(s2.ar,qval)
quant.3a<-quantile(s3.ar,qval)

quant.1b<-quantile(s1.br,qval)
quant.2b<-quantile(s2.br,qval)
quant.3b<-quantile(s3.br,qval)

quantiles<-rbind(quant.1a,quant.2a,quant.3a,quant.1b,quant.2b,quant.3b)
  rownames(quantiles)<-c("s1_dataA","s2_dataA","s3_dataA","s1_dataB","s2_dataB","s3_dataB")
  colnames(quantiles)<-c("0.999","0.990","0.950")

       }   # "if(nci.boo …"   end


# ***************************************************************************** 
# Output
# ***************************************************************************** 
saida<-list(ncol.a,s.obs,quantiles)
names(saida)<-c("NumberOfVariables","s_statistics","Null_Dist_Upper_Quantiles")

return(saida)

}            # end of function






