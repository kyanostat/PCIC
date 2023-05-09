
##Libaries
library(MASS)
library(matrixStats)
library(coda)
library(polynom)
library(ggplot2)

## We use the set-up in Sugiyama et al. (2007)

set.seed(104)
################################################################################
# Parameters
n = 50  # Data sample
m = n
M= 2000 # Posterior sample
p = 1    # dim of X

#sinc function
sinc<-function(x){
  return( sin(pi*x)/(pi*x) )
}

#Data set-up
sigma_Y = 1/4
peak_train  = 0
peak_test   = 0.5
width_train = 1
width_test  = 0.3
X_train = rnorm(n=n,mean=peak_train,sd=width_train) 
Y_train = sinc(X_train) + rnorm(n=n,mean=0,sd=sigma_Y)
X_test  = rnorm(n=m,mean=peak_test,sd=width_test)
Y_test  = sinc(X_test)  + rnorm(n=m,mean=0,sd=sigma_Y)
Xblock_train = cbind(rep(1,n),X_train)
Xblock_test = cbind(rep(1,m),X_test)
################################################################################


################################################################################
#Required functions
densityratio_lambda<-function(x,lambda){
  logp_test =  (dnorm(x,mean=peak_test,sd=width_test,log=TRUE))
  logp_train = (dnorm(x,mean=peak_train,sd=width_train,log=TRUE))
  return( exp( lambda*( logp_test-logp_train) ) )
}
#Log Likelihood
log_p <- function(y,x,beta){
  dif = y-x%*%beta
  ret = log(dnorm(dif,mean=0,sd=sigma_Y))
  return( ret )
}
################################################################################


result <- list() 
PCIC <- list() 
PCIC_body <- list()
PCIC_pen  <- list()
WAIC <- list()
True_Gn <- list()
True_Gn_2 <- list()

counter=0
interval_grid=0.1
interval_upper = 2.

for(lambdas in seq(interval_grid,interval_upper,by=interval_grid)){
  counter=counter+1
  
  ##Set up of Bayes
  #Prior set up
  #Prior for beta N(0, S=tau*I)
  tau = 1
  
  Wblock = diag(sapply(X_train,FUN=densityratio_lambda,lambda=lambdas))
  Sigma_pos = (sigma_Y^2)* ginv( t(Xblock_train) %*% (Wblock) %*% (Xblock_train)   +  (sigma_Y^2 / tau)* diag(2) )
  mu_pos    = Sigma_pos %*% ( t(Xblock_train) %*% (Wblock) %*% Y_train ) / sigma_Y^2
  possample_beta = mvrnorm(n = M, mu=mu_pos, Sigma= Sigma_pos, tol = 1e-06, empirical = FALSE)
  result[[counter]] <- possample_beta
  
  T_n_seq   = seq(1,n)
  V_n_seq = seq(1,n) 
  WAIC_seq_1 = seq(1,n)
  WAIC_seq_2 = seq(1,n)
  for(N in seq(1,n)){
    T_n_seq[N]  = (densityratio_lambda(X_train[N],lambda=1))*(-logSumExp(apply(possample_beta,MARGIN=1,FUN=log_p, 
                                       y = Y_train[N],x = Xblock_train[N,]))-log(1/M))
    
    V_n_seq_1  = mean(apply(possample_beta,MARGIN=1,FUN=log_p,  y = Y_train[N],x = Xblock_train[N,])^2)
    V_n_seq_2  = (mean(apply(possample_beta,MARGIN=1,FUN=log_p,  y = Y_train[N],x = Xblock_train[N,])))^2
    V_n_seq[N] = (densityratio_lambda(X_train[N],lambda=(1+lambdas)))*(V_n_seq_1 - V_n_seq_2)
    
    WAIC_seq_1[N]     = (-logSumExp(apply(possample_beta,MARGIN=1,FUN=log_p,y = Y_train[N],x = Xblock_train[N,]))-log(1/M))
    WAIC_2_1  = mean(apply(possample_beta,MARGIN=1,FUN=log_p,  y = Y_train[N],x = Xblock_train[N,])^2)
    WAIC_2_2  = (mean(apply(possample_beta,MARGIN=1,FUN=log_p,  y = Y_train[N],x = Xblock_train[N,])))^2
    WAIC_seq_2[N]     = (WAIC_2_1 - WAIC_2_2)
    
}
  PCIC[[counter]]      <- mean(T_n_seq + V_n_seq)
  PCIC_body[[counter]] <- mean(T_n_seq)
  PCIC_pen[[counter]]  <- mean(V_n_seq)
  WAIC[[counter]]      <- mean(WAIC_seq_1+WAIC_seq_2)
  
  G_n_seq = seq(1,m)
  G_n_seq_2 = seq(1,m)
  
  for(NN in seq(1,m)){
      
    G_n_seq[NN]   = (-logSumExp(apply(possample_beta,MARGIN=1,FUN=log_p,y = Y_test[NN],x = Xblock_test[NN,]))-log(1/M))
    
    sigma_pred    = sqrt(sigma_Y^2 + c(t(Xblock_test[NN,])%*%Xblock_test[NN,]) * sum(diag(Sigma_pos)) )
    mu_pred       = Xblock_test[NN,]%*%mu_pos
    G_n_seq_2[NN] = -log(dnorm(x=Y_test[NN],mean=mu_pred,sd=sigma_pred))
    
  }
  True_Gn[[counter]]   <- mean(G_n_seq)
  True_Gn_2[[counter]] <- mean(G_n_seq_2)
}

################################################################################

trainingdata_forplot<- data.frame(y=Y_train,x=X_train)
testdata_forplot<- data.frame(y=Y_test,x=X_test)

pcic_forplot<- data.frame(y=unlist(PCIC),    x=seq(interval_grid,interval_upper,by=interval_grid))
waic_forplot<- data.frame(y=unlist(WAIC),    x=seq(interval_grid,interval_upper,by=interval_grid))
gn_forplot<-   data.frame(y=unlist(True_Gn), x=seq(interval_grid,interval_upper,by=interval_grid))

p0 <- ggplot() + theme_bw()
p0 <- p0 + geom_point(data=pcic_forplot,aes(x=x, y=y), stat='identity',alpha = 3/4,colour="blue")
p0 <- p0 + geom_point(data=waic_forplot,aes(x=x, y=y), stat='identity',alpha = 3/4,colour="red")
p0 <- p0 + geom_point(data=gn_forplot,aes(x=x, y=y), stat='identity',alpha = 3/4)
scale_min = min(min(unlist(PCIC)),min(unlist(True_Gn)),min(unlist(WAIC)))
scale_max = max(max(unlist(PCIC)),max(unlist(True_Gn)),max(unlist(WAIC)))
p0 <- p0 + scale_y_continuous(limits = c(-scale_min-0.1,10+0.1))+ xlab(expression(lambda)) + ylab("PCIC, WAIC, and Generalization error")
p0 <- p0 + theme(legend.position = c(1.5, 0.75))
p0 <- p0 + annotate('text', x = 0.95, y = 8., label = "PCIC",color="blue",size=10,family="Times New Roman")
p0 <- p0 + annotate('text', x = 0.7, y = 8., label = "●",color="blue",size=10,family="Times New Roman")
p0 <- p0 + annotate('text', x = 0.95, y = 9., label = "WAIC",size=10,color="red",alpha=3/4,family="Times New Roman")
p0 <- p0 + annotate('text', x = 0.7, y = 9., label = "●",size=10,color="red",alpha=3/4,family="Times New Roman")
p0 <- p0 + annotate('text', x = 1.4, y = 10.0, label = "Generalization error",color="black",alpha=3/4,size=10,family="Times New Roman")
p0 <- p0 + annotate('text', x = 0.7, y = 10.0, label = "●",color="black",alpha=3/4,size=10,family="Times New Roman")
p0 <- p0 + theme(text = element_text(size = 18,family="Times New Roman"))
p0 <- p0 + theme(axis.text = element_text(size = 20))
p0 <- p0 + theme(legend.text = element_text(size = 20))


p00 <- ggplot() + theme_bw()
p00 <- p00 + geom_point(data=pcic_forplot,aes(x=x, y=y), stat='identity',alpha = 3/4,color="blue")
p00 <- p00 + scale_y_continuous(limits = c(0.1,0.5)) +xlab("")+ylab("")
p00 <- p00 + annotate('text', x = 1.35, y = 8., label = "PCIC",color="blue",size=10,family="Times New Roman")
p00 <- p00 + annotate('text', x = 0.8, y = 8., label = "●",color="blue",size=10,family="Times New Roman")
p00 <- p00 + theme(axis.text = element_text(size = 20))
p00 <- p00 + theme(legend.text = element_text(size = 20))


argmin_PCIC <- which.min((unlist(PCIC)))
argmin_WAIC <- which.min((unlist(WAIC)))
p1 <- ggplot() + theme_bw()
p1 <- p1 + geom_point(data=trainingdata_forplot,aes(x=x, y=y), stat='identity',alpha = 3/4)
p1 <- p1 + geom_point(data=testdata_forplot,aes(x=x, y=y), stat='identity',color="chartreuse4",alpha = 3/4)
p1 <- p1 + scale_x_continuous(limits = c(-1, 2))
p1 <- p1 + scale_y_continuous(limits = c(-2, 2))
p1 <- p1 +  xlab(expression(italic(x))) + ylab(expression(italic(y)))
p1 <- p1 + geom_abline(intercept = apply(result[[argmin_PCIC]],MARGIN=2,FUN=mean)[1], slope = apply(result[[argmin_PCIC]],MARGIN=2,FUN=mean)[2],color="blue")
p1 <- p1 + geom_abline(intercept = apply(result[[argmin_PCIC]],MARGIN=2,FUN=quantile,probs=0.025)[1], slope = apply(result[[argmin_PCIC]],MARGIN=2,FUN=quantile,probs=0.05)[2],linetype="dashed",color="blue")
p1 <- p1 + geom_abline(intercept = apply(result[[argmin_PCIC]],MARGIN=2,FUN=quantile,probs=0.975)[1], slope = apply(result[[argmin_PCIC]],MARGIN=2,FUN=quantile,probs=0.95)[2],linetype="dashed",color="blue")
p1 <- p1 + geom_abline(intercept = apply(result[[argmin_WAIC]],MARGIN=2,FUN=mean)[1], slope = apply(result[[argmin_WAIC]],MARGIN=2,FUN=mean)[2],color="red")
p1 <- p1 + geom_abline(intercept = apply(result[[argmin_WAIC]],MARGIN=2,FUN=quantile,probs=0.025)[1], slope = apply(result[[argmin_WAIC]],MARGIN=2,FUN=quantile,probs=0.05)[2],linetype="dashed",color="red")
p1 <- p1 + geom_abline(intercept = apply(result[[argmin_WAIC]],MARGIN=2,FUN=quantile,probs=0.975)[1], slope = apply(result[[argmin_WAIC]],MARGIN=2,FUN=quantile,probs=0.95)[2],linetype="dashed",color="red")
p1 <- p1 + annotate('text', x = 1.5, y = 1.2, label = "Training",size=10,family="Times New Roman")
p1 <- p1 + annotate('text', x = 1.0, y = 1.2, label = "●",size=10,family="Times New Roman")
p1 <- p1 + annotate('text', x = 1.4, y = 1., label = "Test",color="chartreuse4",alpha=4/4,size=10,family="Times New Roman")
p1 <- p1 + annotate('text', x = 1.0, y = 1., label = "●",color="chartreuse4",alpha=4/4,size=10,family="Times New Roman")
p1 <- p1 + theme(text = element_text(size = 18))+ xlab(expression(italic(x))) + ylab(expression(italic(y))) 
p1 <- p1 + theme(axis.text = element_text(size = 20))
p1 <- p1 + theme(legend.text = element_text(size = 20))

print(p0)
