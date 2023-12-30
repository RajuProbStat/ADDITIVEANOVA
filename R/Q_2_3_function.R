library(Iso)
library(matlib)

#' Evaluate the critical point and observed value of the test statistic Q_2 and Q_3
#'
#' @param a number of levels of the row factor A
#' @param b number of levels of the column factor
#' @param size the level of significance
#' @param n matrix having a rows and b columns, and (i,j)-th element is the sample size in the (i,j)-th cell
#' @param mean matrix having a rows and b columns, and (i,j)-th element is the sample mean in the (i,j)-th cell
#' @param var matrix having a rows and b columns, and (i,j)-th element is the sample variance in the (i,j)-th cell
#'
#' @return numeric vector consisting of observed and critical values of Q_2 and Q_3
#' @export
#'
#' @examples
#' Q_2_3(3,2,0.05,rbind(c(15,17),c(16,14),c(18,15)),rbind(c(186.67,332.24),c(261.12,485.50),c(302.44,581.87)),rbind(c(295.10,2374.94),c(372.92,1008.27),c(840.85,661.12)))
Q_2_3<-function(a,b,size,n,mean,var)
{
  B<-10000
  mean_alpha<-rep(NA,a)
  mean_beta<-rep(NA,b)
  mean_boot<-array(NA,dim=c(a,b))
  var_boot<-array(NA,dim=c(a,b))
  mean_alpha_boot<-rep(NA,a)
  mean_beta_boot<-rep(NA,b)
  test1_boot<-rep(NA,B)
  test2_boot<-rep(NA,B)
  for(i in 1:a)
  {
    mean_alpha[i]<-mean(mean[i,])
  }
  for(j in 1:b)
  {
    mean_beta[j]<-mean(mean[,j])
  }
  var1<-rep(0,a-1)
  var2<-rep(0,b-1)
  T1<-rep(NA,a-1)
  T2<-rep(NA,b-1)
  for(i in 1:a-1)
  {
    var1_0<-0
    for(j in 1:b)
    {
      var1[i]<-var1_0+(var[i+1,j]/n[i+1,j])+(var[i,j]/n[i,j])
      var1_0<-var1[i]
    }
    T1[i]<-(mean_alpha[i+1]-mean_alpha[i])/sqrt(var1[i]/b^2)
  }
  for(j in 1:b-1)
  {
    var2_0<-0
    for(i in 1:a)
    {
      var2[j]<-var2_0+(var[i,j+1]/n[i,j+1])+(var[i,j]/n[i,j])
      var2_0<-var2[j]
    }
    T2[j]<-(mean_beta[j+1]-mean_beta[j])/sqrt(var2[j]/a^2)
  }
  test1_ob<-min(c(min(c(T1)),min(c(T2))))      #observed value of the test statistic Q_3
  test2_ob<-min(c(max(c(T1)),max(c(T2))))      #observed value of the test statistic Q_2
  set.seed(17)
  for(l in 1:B)
  {
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        data_boot<-rnorm(n[i,j],0,sqrt(var[i,j]))
        mean_boot[i,j]<-mean(data_boot)
        var_boot[i,j]<-var(data_boot)
      }
      mean_alpha_boot[i]<-mean(mean_boot[i,])
    }
    for(j in 1:b)
    {
      mean_beta_boot[j]<-mean(mean_boot[,j])
    }
    var1_boot<-rep(0,a-1)
    var2_boot<-rep(0,b-1)
    T1_boot<-rep(NA,a-1)
    T2_boot<-rep(NA,b-1)
    for(i in 1:a-1)
    {
      var1_0<-0
      for(j in 1:b)
      {
        var1_boot[i]<-var1_0+(var_boot[i+1,j]/n[i+1,j])+(var_boot[i,j]/n[i,j])
        var1_0<-var1_boot[i]
      }
      T1_boot[i]<-(mean_alpha_boot[i+1]-mean_alpha_boot[i])/sqrt(var1_boot[i]/b^2)
    }
    for(j in 1:b-1)
    {
      var2_0<-0
      for(i in 1:a)
      {
        var2_boot[j]<-var2_0+(var_boot[i,j+1]/n[i,j+1])+(var_boot[i,j]/n[i,j])
        var2_0<-var2_boot[j]
      }
      T2_boot[j]<-(mean_beta_boot[j+1]-mean_beta_boot[j])/sqrt(var2_boot[j]/a^2)
    }
    test1_boot[l]<-min(c(min(c(T1_boot)),min(c(T2_boot))))
    test2_boot[l]<-min(c(max(c(T1_boot)),max(c(T2_boot))))
  }
  critical_test1<-quantile(test1_boot,1-size)                 #critical point of the test statistic Q_3
  critical_test2<-quantile(test2_boot,1-size)                 #critical point of the test statistic Q_2
  cat("critical point of Q_2, observed value of Q_2, critical point of Q_3 and observed value of Q_3 are\n")
  return(c(critical_test2,test2_ob,critical_test1,test1_ob))
}
