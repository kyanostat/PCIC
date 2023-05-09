##Libaries
library(MASS)
library(matrixStats)
library(coda)
library(polynom)
library(ggplot2)
## Comparison between PCIC and Cp
## We use the set-up in Baba, Kanamori, and Ninomiya Biometrika (2017)

set.seed(122)


################################################################################
# Parameters
n = 50  # Data sample
M= 1000 # Posterior sample
p = 2    # dim of X
H = 6    # num of assingments
b0 = c(1, 1, 0.5) # True of b0

##Set up of Bayes
#Prior set up
#Prior for beta N(0, S/tau)
assumed_p = p+1
S = 0.001*diag(c(rep(1,assumed_p),rep(0,(p+1-assumed_p))))

################################################################################
#Output
#Cp
#WCPi #Sample-wise WCp
#UCPi #Sample-wise UCp
#pen_UCPi #Sample-wise penalty in UCp
#body_UCPi #Sample-wise main term in UCp
#pen_WCPi #Sample-wise penalty in WCp
#body_WCPi #Sample-wise main term in WCp
#WWAIC
#T_n_seq  # Tn for weighted kernel
#UT_n_seq  # Tn for unweighted kernel
#V_n_seq  # Vn for weighted kernel
#UV_n_seq  # Vn for unweighted kernel
#WWWAIC = T_n_seq + V_n_seq
#UWWAIC = UT_n_seq + UV_n_seq
#LOOCV
#WLOOCV_n_seq  # CV for weight
#ULOOCV_n_seq  # CV for unweight

################################################################################
#######
##functions required
#Generate X
polygen = function(h,p){
  ret = h^0
  for(i in c(1:p)){
    ret = c(ret,h^i)
  }
  return(ret)
}

Xgen = function(H,p){
  ret = polygen(1,p)
  for(h in c(2:H)){
    ret = rbind(ret,polygen(h,p))
  }
  return(ret)
}
#Create propensity score
propensityscore = function(z){
  ret = exp(alpha*z)
  return(ret/sum(ret))
}
#Create assignments
createTis = function(prob){
  ret = rmultinom(1,size=1,prob=prob)
  return(ret)
}


#Log Likelihood for the training
log_s_train <- function(theta,y,x,w){
  #theta = (beta,tau)
  dim = length(x)
  ret1 = (-rev(theta)[1]/2) * (y-x%*%theta[c(1:dim)])^2
  ret2 = (-1/2)*log(2*pi/rev(theta)[1])
  return( w*(ret1 + ret2) )
}

#Log Likelihood for the test
log_h_test <- function(theta,y,x){
  #theta = (beta,tau)
  dim = length(x)
  ret1 = (-rev(theta)[1]/2) * (y-x%*%theta[c(1:dim)])^2
  ret2 = (-1/2)*log(2*pi/rev(theta)[1])
  return( (ret1 + ret2) )
}

################################################################################

################################################################################
# Generate data
X_nonortho = Xgen(H,p) #Covariate  X
Z = runif(n,-sqrt(3),sqrt(3)) # Confounder Z
epsilon = rnorm(n,0,1) # noise independent of Z
varepsilon = Z+epsilon # noise in Y
Ys = c(X_nonortho%*%b0)%o%rep(1,n)  + (c(1,1,1,1,1,1) %o% varepsilon) # Y vector
alpha = c(0, 0.8, 1.0, 0.9, 0.7, 0.6) # param of propensity score
eis = apply(array(Z),MARGIN=1,FUN = propensityscore) # generate propensity score
Tis = apply(eis,MARGIN=2, FUN=createTis) # generate assingments
Wis = Tis / eis # generate W (weight t/e)


##Block diagonalization of Y, X, W
Yblock = c(t(Ys)) #Create block Y
r_box = list()
for(h in c(1:H)){
  r_box[[h]] = diag(Wis[h,])
  if(h==1){
    Wblock = rbind(r_box[[h]],0*matrix(0,nrow=n*(H-1),ncol=n))
  }else{
    intermid=rbind(0*matrix(0,nrow=n*(h-1),ncol=n),r_box[[h]],0*matrix(0,nrow=n*(H-h),ncol=n))
    Wblock = cbind(Wblock,intermid)
  }  
}
Xblock_pre = 0*diag((p+1))
for(h in c(1:H)){
  Xblock_pre = Xblock_pre + (X_nonortho [h,]) %o% X_nonortho [h,]
}
A = chol(Xblock_pre)
X = X_nonortho%*%ginv(A)
Xblock = matrix((rep(1,n) %o% X),ncol=(p+1),nrow=n*H)
################################################################################
###############################################################################
#Generate posterior samples
Sigma_pos = ginv( t(Xblock)%*%Wblock%*%Xblock  + S)
mu_pos    = Sigma_pos %*% ( t(Xblock) %*% Wblock %*% Yblock )

possample_beta = mvrnorm(n = M, mu=mu_pos, Sigma= (var(varepsilon))*Sigma_pos, tol = 1e-06, empirical = FALSE)
possample_tau  = rep(1/var(varepsilon) ,M)
possample_theta = cbind(possample_beta,possample_tau)
###############################################################################
###############################################################################
##Observation
Y_obs = apply(Tis*Ys,MARGIN=2,FUN=sum)
X_obs = t(Tis) %*% X
###############################################################################

T_n_seq = seq(1,n) # Tn for weighted kernel
UT_n_seq = seq(1,n) # Tn for unweighted kernel
V_n_seq = seq(1,n) # Vn for weighted kernel
UV_n_seq = seq(1,n) # Vn for unweighted kernel


for(N in c(1:n)){
  T_n_seq[N] = (sum(Wis[,N]))*(-logSumExp(apply(possample_theta,MARGIN=1,FUN=log_h_test, 
                                y = Y_obs[N],x = X_obs[N,]))-log(1/M))
  UT_n_seq[N] = (1)*(-logSumExp(apply(possample_theta,MARGIN=1,FUN=log_h_test, 
                                 y = Y_obs[N],x = X_obs[N,]))-log(1/M))
  
  V_n_seq_1 = mean(apply(possample_theta,MARGIN=1,FUN=log_s_train,  y = Y_obs[N],x = X_obs[N,], w=sum(Wis[,N]))^2)
  V_n_seq_2 = (mean(apply(possample_theta,MARGIN=1,FUN=log_s_train,  y = Y_obs[N],x = X_obs[N,], w=sum(Wis[,N]))))^2
  V_n_seq[N] = V_n_seq_1 - V_n_seq_2
  
  UV_n_seq_1 = mean(apply(possample_theta,MARGIN=1,FUN=log_s_train,  y = Y_obs[N],x = X_obs[N,], w=1)^2)
  UV_n_seq_2 = (mean(apply(possample_theta,MARGIN=1,FUN=log_s_train,  y = Y_obs[N],x = X_obs[N,], w=1)))^2
  UV_n_seq[N] = sum(Wis[,N])* ( UV_n_seq_1 - UV_n_seq_2 )
}

WPCIC = T_n_seq + V_n_seq
UPCIC = UT_n_seq + UV_n_seq


################################################################################
#Cp (UCp and WCp in Baba, Kanamori, and Ninomiya Biometrika 2017)

##dim calculation
dimi_WCP = function(y,x,w,ipw){
  ret = (w^2)* (t(x)%*%(x)) * (y-x%*% ipw)^2
  return(ret)
}

mu_IPW    = ginv( t(Xblock)%*%Wblock%*%Xblock  ) %*% ( t(Xblock) %*% Wblock %*% Yblock )
WCPi = seq(1,n)
UCPi = seq(1,n)
pen_UCPi = seq(1,n)
body_UCPi = seq(1,n)
pen_WCPi = seq(1,n)
body_WCPi = seq(1,n)
for(N in c(1:n)){
  body_UCPi[N]  =  (1/(2*var(varepsilon)))* (Y_obs[N] - X_obs[N,] %*% mu_IPW) * (Y_obs[N] - X_obs[N,] %*% mu_IPW)  + (1/2)*log(2*var(varepsilon)*pi)
  pen_UCPi[N]   =  (p+1)/n
  body_WCPi[N] = (1/(2*var(varepsilon)))* (Y_obs[N] - X_obs[N,] %*% mu_IPW) * sum(Wis[,N]) *  (Y_obs[N] - X_obs[N,] %*% mu_IPW) + (sum(Wis[,N])/2)*log(2*var(varepsilon)*pi)
  pen_WCPi[N] = dimi_WCP(y=Y_obs[N],x=X_obs[N,],w=sum(Wis[,N]),ipw=mu_IPW)/(var(varepsilon)* n )
  WCPi[N] = body_WCPi[N] + pen_WCPi[N]
  UCPi[N] = body_UCPi[N] + pen_UCPi[N]
}
################################################################################


W_forplot<- data.frame(y=WPCIC,x=WCPi)
Wpen_forplot<- data.frame(y=V_n_seq,x=pen_WCPi)

p1 <- ggplot() + theme_bw()
p1 <- p1 + geom_point(data=W_forplot,aes(x=x, y=y), stat='identity',alpha = 3/4,color="blue")
scale_min = 0
scale_max = 20
p1 <- p1 + scale_y_continuous(limits = c(scale_min,scale_max)) 
p1 <- p1 + scale_x_continuous(limits = c(scale_min,scale_max)) 
p1 <- p1 + xlab("PCICw") + ylab("wCp")
p1 <- p1 + geom_abline(intercept = 0, slope = 1,alpha=0.5)
p1 <- p1 + theme(text = element_text(size = 20,family="Times New Roman"))
p1 <- p1 + theme(axis.text = element_text(size = 20))
p1 <- p1 + theme(legend.text = element_text(size = 20))

p2 <- ggplot() + theme_bw()
p2 <- p2 + geom_point(data=Wpen_forplot,aes(x=x, y=y), stat='identity',alpha = 3/4,color="blue")
scale_min = 0
scale_max = 2.2
p2 <- p2 + scale_y_continuous(limits = c(scale_min,scale_max)) 
p2 <- p2 + scale_x_continuous(limits = c(scale_min,scale_max)) 
p2 <- p2 + xlab("Bias correction term of PCICw") + ylab("Bias correction term of wCp")
p2 <- p2 + geom_abline(intercept = 0, slope = 1,alpha=0.5)
p2 <- p2 + theme(text = element_text(size = 20,family="Times New Roman"))
p2 <- p2 + theme(axis.text = element_text(size = 20))
p2 <- p2 + theme(legend.text = element_text(size = 20))

print(p1)

