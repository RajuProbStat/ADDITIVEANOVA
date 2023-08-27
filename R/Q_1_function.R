library(Iso)
library(matlib)

#' Evaluate the critical point and observed value of the test statistic Q_1
#'
#' @param a number of levels of the row factor A
#' @param b number of levels of the column factor B
#' @param size the level of significance
#' @param n matrix having a rows and b columns, and (i,j)-th element is the sample size in the (i,j)-th cell
#' @param mean matrix having a rows and b columns, and (i,j)-th element is the sample mean in the (i,j)-th cell
#' @param var matrix having a rows and b columns, and (i,j)-th element is the sample variance in the (i,j)-th cell
#'
#' @return numeric vector consisting of observed and critical values of Q_1
#' @export
#'
#' @examples
#' Q_1(3,2,0.05,rbind(c(15,17),c(16,14),c(18,15)),rbind(c(186.67,332.24),c(261.12,485.50),c(302.44,581.87)),rbind(c(295.10,2374.94),c(372.92,1008.27),c(840.85,661.12)))
Q_1<-function(a,b,size,n,mean,var){
  func_MLE<-function(a,b,n,mean,var){
    ########################################
    #MLEs of the parameters under \Omega_{0A}
    mean_alpha_0A<-rep(NA,a)
    mean_beta_0A<-rep(NA,b)
    for(i in 1:a)
    {
      mean_alpha_0A[i]<-mean(mean[i,])
    }
    for(j in 1:b)
    {
      mean_beta_0A[j]<-mean(mean[,j])
    }
    beta_0_0A<-mean_beta_0A                      #initialization of beta
    var_0_0A<-array(NA,c(a,b))
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        var_0_0A[i,j]<-var[i,j]                 #initialization of cell variances
      }
    }
    beta_1_0A<-rep(NA,b)
    var_1_0A<-array(NA,c(a,b))
    repeat
    {
      for(j in 1:b)
      {
        sum_1_A<-0
        sum_2_A<-0
        for(i in 1:a)
        {
          sum_1_A<-sum_1_A+(n[i,j]/var_0_0A[i,j])
          sum_2_A<-sum_2_A+(n[i,j]/var_0_0A[i,j])*mean[i,j]
        }
        beta_1_0A[j]<-sum_2_A/sum_1_A                           #updated value of beta
      }
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_1_0A[i,j]<-var[i,j]+(mean[i,j]-beta_1_0A[j])^2    #updated values of cell variances
        }
      }
      if(max(abs(beta_1_0A-beta_0_0A))<0.00001 & norm(var_1_0A-var_0_0A, type="M")<0.00001) #stopping criteria
      {
        break
      }
      beta_0_0A<-beta_1_0A
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_0_0A[i,j]<-var_1_0A[i,j]
        }
      }
    }
    ########################################
    ########################################
    #MLEs of the parameters under \Omega_{A}
    alpha_0_A<-rep(NA,a)
    beta_0_A<-rep(NA,b)
    var_0_A<-array(NA,c(a,b))
    alpha_1_A<-rep(NA,a)
    beta_1_A<-rep(NA,b)
    var_1_A<-array(NA,c(a,b))
    for(i in 1:a)
    {
      alpha_0_A[i]<-mean(mean[i,])             #initialization of alpha
    }
    for(j in 1:b)
    {
      beta_0_A[j]<-mean(mean[,j])-mean(mean)   #initialization of beta
    }
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        var_0_A[i,j]<-var[i,j]                #initialization of cell variances
      }
    }
    repeat
    {
      u_A<-array(NA,c(a,b))
      w_A<-rep(NA,a)
      x_A<-rep(NA,a)
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          u_A[i,j]<-(n[i,j]/var_0_A[i,j])
        }
      }

      for(i in 1:a)
      {
        sum_3_A<-0
        for(j in 1:b)
        {
          sum_3_A<-sum_3_A+u_A[i,j]
        }
        w_A[i]<-sum_3_A
        sum_4_A<-0
        for(j in 1:b)
        {
          sum_4_A<-sum_4_A+u_A[i,j]*mean[i,j]
        }
        sum_5_A<-0
        for(j in 1:(b-1))
        {
          sum_5_A<- sum_5_A+(u_A[i,b]-u_A[i,j])*beta_0_A[j]
        }
        x_A[i]<-(1/w_A[i])*(sum_4_A+sum_5_A)
      }
      alpha_1_A<-pava(x_A, w_A, decreasing=FALSE, long.out=FALSE, stepfun=FALSE) #updated value of alpha
      t_A<-rep(NA,b-1)
      for(j in 1:(b-1))
      {
        sum_6_A<-0
        for(i in 1:a)
        {
          sum_6_A<-sum_6_A+u_A[i,j]*(mean[i,j]-alpha_1_A[i])-u_A[i,b]*(mean[i,b]-alpha_1_A[i])
        }
        t_A[j]<-sum_6_A
      }
      Q_1A<-rep(NA,b-1)
      for(j in 1:(b-1))
      {
        Q_1A[j]<-a*mean(u_A[,j])
      }
      Q_A<-diag(Q_1A,nrow=(b-1),ncol=(b-1))+a*mean(u_A[,b])*matrix(1,nrow=b-1,ncol=b-1)
      beta_1_A_old<-Ginv(Q_A)%*%t_A
      for(j in 1:(b-1))
      {
        beta_1_A[j]<-beta_1_A_old[j]                                             #updated value of beta
      }
      beta_1_A[b]<- -sum(beta_1_A_old)
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_1_A[i,j]<-var[i,j]+(mean[i,j]-alpha_1_A[i]-beta_1_A[j])^2          #updated values of cell variances
        }
      }
      if(max(abs(alpha_1_A-alpha_0_A))<0.00001 & max(abs(beta_1_A-beta_0_A))<0.00001) #stopping criteria
      {
        break
      }
      alpha_0_A<-alpha_1_A
      beta_0_A<-beta_1_A
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_0_A[i,j]<-var_1_A[i,j]
        }
      }
    }
    #######################################
    #######################################
    #MLEs of the parameters under \Omega_0B
    mean_alpha_0B<-rep(NA,a)
    mean_beta_0B<-rep(NA,b)
    for (i in 1:a)
    {
      mean_alpha_0B[i]<-mean(mean[i,])       #initialization of alpha
    }
    for (j in 1:b)
    {
      mean_beta_0B[j]<-mean(mean[,j])
    }
    alpha_0_0B<-mean_alpha_0B
    var_0_0B<-array(NA,c(a,b))
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        var_0_0B[i,j]<-var[i,j]                                  #initialization of cell variances
      }
    }
    alpha_1_0B<-rep(NA,a)
    var_1_0B<-array(NA,c(a,b))
    repeat
    {
      for(i in 1:a)
      {
        sum_1_B<-0
        sum_2_B<-0
        for(j in 1:b)
        {
          sum_1_B<-sum_1_B+(n[i,j]/var_0_0B[i,j])
          sum_2_B<-sum_2_B+(n[i,j]/var_0_0B[i,j])*mean[i,j]
        }
        alpha_1_0B[i]<-sum_2_B/sum_1_B                           #updated value of alpha
      }
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_1_0B[i,j]<-var[i,j]+(mean[i,j]-alpha_1_0B[i])^2    #updated values of cell variances
        }
      }
      if(max(abs(alpha_1_0B-alpha_0_0B))<0.00001 & norm(var_1_0B-var_0_0B, type="M")<0.00001) #stopping criteria
      {
        break
      }
      alpha_0_0B<-alpha_1_0B
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_0_0B[i,j]<-var_1_0B[i,j]
        }
      }
    }
    ######################################
    ######################################
    #MLEs of the parameters under Omega_B
    alpha_0_B<-rep(NA,a)
    beta_0_B<-rep(NA,b)
    var_0_B<-array(NA,c(a,b))
    alpha_1_B<-rep(NA,a)
    beta_1_B<-rep(NA,b)
    var_1_B<-array(NA,c(a,b))
    for(j in 1:b)
    {
      beta_0_B[j]<-mean(mean[,j])                   #initialization of beta
    }
    for(i in 1:a)
    {
      alpha_0_B[i]<-mean(mean[i,])-mean(mean)       #initialization of alpha
    }
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        var_0_B[i,j]<-var[i,j]                      #initialization of cell variances
      }
    }
    repeat
    {
      u_B<-array(NA,c(a,b))
      w_B<-rep(NA,b)
      x_B<-rep(NA,b)
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          u_B[i,j]<-(n[i,j]/var_0_B[i,j])
        }
      }
      for(j in 1:b)
      {
        sum_3_B<-0
        for(i in 1:a)
        {
          sum_3_B<-sum_3_B+u_B[i,j]
        }
        w_B[j]<-sum_3_B
        sum_4_B<-0
        for(i in 1:a)
        {
          sum_4_B<-sum_4_B+u_B[i,j]*mean[i,j]
        }
        sum_5_B<-0
        for(i in 1:(a-1))
        {
          sum_5_B<- sum_5_B+(u_B[a,j]- u_B[i,j])*alpha_0_B[i]
        }
        x_B[j]<-(1/w_B[j])*(sum_4_B+sum_5_B)
      }
      beta_1_B<-pava(x_B, w_B, decreasing=FALSE, long.out=FALSE, stepfun=FALSE)  #updated value of beta
      t_B<-rep(NA,a-1)
      for(i in 1:(a-1))
      {
        sum_6_B<-0
        for(j in 1:b)
        {
          sum_6_B<-sum_6_B+u_B[i,j]*(mean[i,j]-beta_1_B[j])-u_B[a,j]*(mean[a,j]-beta_1_B[j])
        }
        t_B[i]<-sum_6_B
      }
      Q_1B<-rep(NA,a-1)
      for(i in 1:(a-1))
      {
        Q_1B[i]<-b*mean(u_B[i,])
      }
      Q_B<-diag(Q_1B,nrow=(a-1),ncol=(a-1))+b*mean(u_B[a,])*matrix(1,nrow=a-1,ncol=a-1)
      alpha_1_B_old<-Ginv(Q_B)%*%t_B
      for(i in 1:(a-1))
      {
        alpha_1_B[i]<-alpha_1_B_old[i]                                           #updated value of alpha
      }
      alpha_1_B[a]<- -sum(alpha_1_B_old)
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_1_B[i,j]<-var[i,j]+(mean[i,j]-beta_1_B[j]-alpha_1_B[i])^2          #updated values of cell variance
        }
      }
      if(max(abs(alpha_1_B-alpha_0_B))<0.00001 & max(abs(beta_1_B-beta_0_B))<0.00001) #stopping criteria
      {
        break
      }
      alpha_0_B<-alpha_1_B
      beta_0_B<-beta_1_B
      for(i in 1:a)
      {
        for(j in 1:b)
        {
          var_0_B[i,j]<-var_1_B[i,j]
        }
      }
    }
    ######################################
    prod_A<-1
    prod_B<-1
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        prod_A<-prod_A*(var_1_A[i,j]/var_1_0A[i,j])^(n[i,j]/2)
        prod_B<-prod_B*(var_1_B[i,j]/var_1_0B[i,j])^(n[i,j]/2)
      }
    }
    return(max(prod_A,prod_B))
  }
  B<-10000
  ob_value<-func_MLE(a,b,n,mean,var)               #observed value of the test statistic Q_1
  mean_boot<-array(NA,c(a,b))
  var_boot<-array(NA,c(a,b))
  MLRT_boot<-rep(NA,B)
  for(k in 1:B)
  {
    for(i in 1:a)
    {
      for (j in 1:b)
      {
        data_boot<-rnorm(n[i,j],0,(var[i,j])^(1/2))
        mean_boot[i,j]<-mean(data_boot)
        var_boot[i,j]<-var(data_boot)
      }
    }
    MLRT_boot[k]<-func_MLE(a,b,n,mean_boot,var_boot)
  }
  # plot(MLRT_boot)
  critical_value<-quantile(MLRT_boot,size)       #critical value of the test statistic Q_1
  cat("observed value of Q_1, critical point of Q_1 are\n")
  return(c(ob_value,critical_value))
}



