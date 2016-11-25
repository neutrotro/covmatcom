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

#' These s "squared" statistics are sensitive to the distribution of the        
#' differentiation across eigenvectors. The maximum s1 value (= 1) occurs in   
#' the situation of "maximum concentration" of the differentiation, where one  
#' eigenvector explains all variation in the first dataset and another,        
#' orthogonal eigenvector explains all variation in the second dataset.         
#' Situations with more than one greater than zero eigenvalues per dataset     
#' result in s1 values <1, even when the two eigenvector sets are mutually     
#' orthogonal. This is due to the use of squared values in the calculation of  
#' the s set of statistics. I am grateful to Dr Evan Remington for calling my  
#' attention to this algorithm limitation.                                     
#'                                                                             
#' As an alternative for cases in which detecting matrix orthogonality is      
#' important, covmatcom calculates absolute value-based versions of the same     
#' statistics (s1Abs, s2Abs, s3Abs) not dependent on eigenvalue distribution  
#' and producing maximum s1Abs values (s1Abs = 1) in all orthogonal            
#' situations. The s1Abs values are decomposed in matrix orientation (s2Abs)   
#' and matrix shape contributions (s3Abs) as follows: 
#' 
#' s1Abs = sum |vi11 - vi21| + sum |vi12 - vi22|
#'
#' s2Aux = sum |vi11 - vi22| + sum |vi12 - vi21|
#' s3Aux = sum |vi11 - vi12| + sum |vi21 - vi22|
#'
#' s2Abs = s1Abs s2Aux / (s2Aux + s3Aux)
#' s3Abs = s1Abs s3Aux / (s2Aux + s3Aux)
#'
#' where vijk is the amount of variance of eigenvector i from dataset k when 
#' applied to data in dataset j, and sum is a summation from i=1 to i=n, the
#' number of variables and eigenvectors in each compared dataset.
#'_____________________________________________________________________________
#'
#' @param data.a Dataset to calculate the first cov matrix to be compared
#' @param data.b Dataset to calculate the second matrix to be compared
#' @param nci.boo number of bootstrap resamples to obtain upper confidence
#'        intervals. If = 0, no C.I. are calculated. Defaults to  = 0.
#' @keywords covMatrix, comparison, orientation, shape
#' @author Carlos Garcia, \email{carlos.garcia.suarez@usc.es}
#' @references Garcia C. 2012. A simple procedure for the comparison of  
#'    covariance matrices. BMC Evolutionary Biology 12:1-17
#'
#' @export
#' @return A list including "NumberOfVariables", "s_squared_statistics", 
#'   "s_squared_Null_Dist_Upper_Quantiles", "s_absolute_statistics" and
#'   "s_absolute_Null_Dist_Upper_Quantiles")

#'    NOTE: The function obtains these null distribution quantiles by 
#'  calculating the s and s_Abs statistics for pairs of bootstrap resamples
#'  taken from a single data set. A set of s and s_Abs statistics' quantiles
#'  is obtained for dataA and another for dataB.



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
  # - "NumberOfVariables": Number of variables 
  # - "s_squared_statistics": s1, s2 and s3 statistics values 
  # - "s_squared_Null_Dist_Upper_Quantiles": 0.999, 0.990 and 0.950 quantiles
  # -   for the three statistics from null distributions obtained by 
  # -   bootstrap resampling within each data set
  # - "s_absolute_statistics": s1Abs, s2Abs and s3Abs statistics values
  # - "s_absolute_Null_Dist_Upper_Quantiles": 0.999, 0.990 and 0.950 quantiles
  # -   for the three statistics from null distributions obtained by 
  # -   bootstrap resampling within each data set


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


# **** "squared" statistics
s1o<-2*(sum(dif.a^2)+sum(dif.b^2)) 
s2o<-sum(diag^2) 
s3o<-sum(hori^2) 

# All "squared statistics" calculations include divisions by 8 to obtain 
# 0 to 1 values for s1. 
# This is because the maximum possible values for sum(dif.a²) and sum(dif.b²)
# are = 2 (explained variances vexp. are expressed as proportions). These
# values would occur if all variance in each sample were explained by single,
# different eigenvectors, so that some, say dif.a[ip] would be equal to +1
# and another, dif.a[ip'], to -1. Thus the maximum possible value for 
# sum(dif.a²) is = 2 (and also =2 for sum(dif.b²)), and is = 2*(2+2) = 8 for
# s1o.

s1<-s1o/8 # total differentiation "squared" statistic
s2<-s2o/8 # orientation contribution "squared" statistic
s3<-s3o/8 # shape contribution "squared" statistic

s.obs<-c(s1,s2,s3)
names(s.obs)<-c("s1-total","s2-orientation","s3-shape")


# **** "absolute" statistics
s1Abs<-(sum(abs(dif.a))+sum(abs(dif.b)))/4

# All s1Abs calculations include divisions by 4 to obtain 0 to 1 
# values for s1Abs. 
# This is because the maximum possible values for sum(abs(dif.a)) and for
# sum(abs(dif.b)) are = 2, i.e., when all variance was explained by different
# eigenvectors in each dataset.

s2oAux<-sum(abs(diag)) 
s3oAux<-sum(abs(hori)) 

TotNoDif<-s2oAux+s3oAux

if(TotNoDif ==0) {s2Abs<-0;s3Abs<-0} else {
s2Abs<-s1Abs*(s2oAux /(TotNoDif))
s3Abs<-s1Abs*(s3oAux /(TotNoDif))}
# Note that s2Abs and s3Abs need no final corrections because they 
# were calculated as proportions

sAbs.obs<-c(s1Abs,s2Abs,s3Abs)
names(sAbs.obs)<-c("s1Abs-total","sAbs2-orientation","sAbs3-shape")



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

s1Abs.ar<-numeric(nci.boo ) 
s2Abs.ar<-numeric(nci.boo )  
s3Abs.ar<-numeric(nci.boo ) 
s1Abs.br<-numeric(nci.boo )  
s2Abs.br<-numeric(nci.boo )  
s3Abs.br<-numeric(nci.boo ) 


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

# All "squared statistics" calculations include divisions by 8 to obtain 
# 0 to 1 values for s1.apr and s1.bpr. 
# This is because the maximum possible values for sum(dif.a1r²),
# sum(dif.a2r²), sum(dif.b1r²) and sum(dif.b2r²) are = 2 (explained variances
# vexp. are expressed as proportions). These values would occur if all
# variance in each sample were explained by single different eigenvectors, so
# that some, say dif.a1r[ip] would be equal to +1 and another, dif.a1r[ip'],
# to -1. Thus the maximum possible value for sum(dif.a1r²) is = 2 (and also =2
# for sum(dif.a21r²)), and is = 2*(2+2) = 8 for s1apr and s1bpr.

s1.ar[iboo]<-s1.apr/8 # total differentiation "squared" statistic
s2.ar[iboo]<-s2.apr/8 # orientation contribution "squared" statistic
s3.ar[iboo]<-s3.apr/8 # shape contribution "squared" statistic

s1.br[iboo]<-s1.bpr/8 # total differentiation "squared" statistic
s2.br[iboo]<-s2.bpr/8 # orientation contribution "squared" statistic
s3.br[iboo]<-s3.bpr/8 # shape contribution "squared" statistic



# **** "absolute" statistics

# All s1Abs.apr and s1Abs.bpr calculations include divisions by 4 to obtain 0
# to 1 values for these statistics. 
# This is because the maximum possible values for sum(abs(dif.a1r)),
# sum(abs(dif.a2r)), sum(abs(dif.b1r)) and sum(abs(dif.b2r)) are = 2, i.e.,
# when all variance was explained by different eigenvectors in each dataset.

s1Abs.apr<-(sum(abs(dif.a1r))+sum(abs(dif.a2r)))/4
s2oAux.apr<-sum(abs(diag.ar)) 
s3oAux.apr<-sum(abs(hori.ar)) 

s1Abs.bpr<-(sum(abs(dif.b1r))+sum(abs(dif.b2r)))/4
s2oAux.bpr<-sum(abs(diag.br)) 
s3oAux.bpr<-sum(abs(hori.br)) 




TotNoDif.apr<-s2oAux.apr+s3oAux.apr
TotNoDif.bpr<-s2oAux.bpr+s3oAux.bpr

if(TotNoDif.apr ==0) {s2oAbs.apr<-0;s3oAbs.apr<-0} else {
s2Abs.apr<-s1Abs.apr*(s2oAux.apr /(TotNoDif.apr))
s3Abs.apr<-s1Abs.apr*(s3oAux.apr /(TotNoDif.apr))}
 
if(TotNoDif.bpr ==0) {s2Abs.bpr<-0;s3Abs.bpr<-0} else {
s2Abs.bpr<-s1Abs.bpr*(s2oAux.bpr /(TotNoDif.bpr))
s3Abs.bpr<-s1Abs.bpr*(s3oAux.bpr /(TotNoDif.bpr))}

# Note that s2Abs.apr,s3Abs.apr s2Abs.bpr and s3Abs.bpr need no further
# corrections because they were calculated as proportions

s1Abs.ar[iboo]<-s1Abs.apr # total differentiation "absolute" statistic
s2Abs.ar[iboo]<-s2Abs.apr # orientation contribution "absolute" statistic
s3Abs.ar[iboo]<-s3Abs.apr # shape contribution "absolute" statistic

s1Abs.br[iboo]<-s1Abs.bpr # total differentiation "absolute" statistic
s2Abs.br[iboo]<-s2Abs.bpr # orientation contribution "absolute" statistic
s3Abs.br[iboo]<-s3Abs.bpr # shape contribution "absolute" statistic

  
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

s.quantiles<-rbind(quant.1a,quant.2a,quant.3a,quant.1b,quant.2b,quant.3b)

  rownames(s.quantiles)<-c("s1_dataA","s2_dataA","s3_dataA","s1_dataB",
           "s2_dataB","s3_dataB")

  colnames(s.quantiles)<-c("0.999","0.990","0.950")

quantAbs.1a<-quantile(s1Abs.ar,qval)
quantAbs.2a<-quantile(s2Abs.ar,qval)
quantAbs.3a<-quantile(s3Abs.ar,qval)

quantAbs.1b<-quantile(s1Abs.br,qval)
quantAbs.2b<-quantile(s2Abs.br,qval)
quantAbs.3b<-quantile(s3Abs.br,qval)

sAbs.quantiles<-rbind(quantAbs.1a,quantAbs.2a,quantAbs.3a,quantAbs.1b,
  quantAbs.2b,quantAbs.3b)

  rownames(sAbs.quantiles)<-c("s1Abs_dataA","s2Abs_dataA","s3Abs_dataA",
    "s1Abs_dataB","s2Abs_dataB","s3Abs_dataB")

  colnames(sAbs.quantiles)<-c("0.999","0.990","0.950")


       }   # "if(nci.boo …"   end

if(nci.boo == 0) {s.quantiles<-NULL;sAbs.quantiles<-NULL}
  

# ***************************************************************************** 
# Output
# ***************************************************************************** 
saida<-list(ncol.a,s.obs,s.quantiles,sAbs.obs,sAbs.quantiles)
names(saida)<-c("NumberOfVariables","s_squared_statistics",
  "s_squared_Null_Dist_Upper_Quantiles","s_absolute_statistics",
  "s_absolute_Null_Dist_Upper_Quantiles")

return(saida)

}            # end of function






