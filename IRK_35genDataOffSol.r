###################################################
# Generate data for wombling starting away from the solution

IRK_35genDataOffSol<- function(I = 5, no.female = 3, s2s = 1, R=2,
                         s2 = 0.5, b0 = 3, b1=0.5, rho=0.9, mcmc=10000){
  K<- 35
  
  # check
  no.obs<- I*K*R
  no.obs
  
  # how will these people be allocated gender status?
  x1 = c(rep(1, no.female*K*R), rep(0, (I-no.female)*R*K) ) # covariate values
  #x1.short = c(rep(1, no.female*R), rep(0, (I-no.female)*K))
  x1.short = c(rep(1, no.female), rep(0, I-no.female))
  
  load("contig_35Reg.Rdata")
  W = contig_35Reg
  #load("strong_spat35Reg.Rdata")
  #W = strong_spat35Reg
  diag(W)<- 0  
  W
  
  # generate data for global mean b0 intercept wombled model
  # with 0 mean random effects
  
  load("RandomW_35region.Rdata")
  start.W = RandomW_35region[,,2]  # try second configuration - looks like random dots
                                   # to recover strong spat matrix

  #start.W = matrix(w, K,K)
  Q = rho*(diag(K)*rowSums(W) - W) + (1-rho)*diag(K)
  
  spatDat<- mvrnorm(I, mu = rep(0, K), Sigma = s2s*solve(Q) )
  library(mefa)
  vec.d<- as.matrix(rep(data.frame(spatDat), each = R))# repeat every random effect R times
  
  vec.data<- as.vector(t(vec.d))
  
  y = b0 + b1*x1 + vec.data + rnorm(I*K*R, mean = 0, sd = sqrt(s2))
  y.vec = y
  
  ##
  dat<- data.frame(y= y.vec, region = factor(rep(1:K, I*R)), person = factor(rep(1:I, each = R*K)))
  
  # windows()
  p.dat<- ggplot(dat, aes(x=region, y=y, colour = region)) + geom_point()+
    ggtitle(paste("CAR simulated data for ", I, "ppl with ", R, "replicates")) + theme_bw() + facet_wrap(~person, ncol = 5)
  
  tR<- matrix(y, nrow = K, ncol = I*R)
  y.mat<- t(tR)
  
  
  Wkeep<-array(0, dim = c(K,K, mcmc)); Wkeep[,,1]<- start.W;   
  beta0 <- beta1<- sigma2<- sigma2.s<- rep(0, mcmc) 
  
  beta0[1] = b0; beta1[1] = b1; 
  sigma2[1] = s2; sigma2.s[1] = s2s; 
  beta.vec = c(b0,b1)
  
  # need to have data in both matrix and vector form
  X.mat = matrix(0, nrow = I, ncol = 2)
  X.mat[, 1]<- rep(1, I)
  X.mat[, 2]<- x1.short
  
  X.short = matrix(0, nrow = I*R, ncol = 2) 
  X.short[, 1]<- rep(1, I*R)
  X.short[, 2]<- x1.short
  #  X.vec<- matrix(0, nrow = I*K*R, ncol = length(beta.vec))
  #  X.vec[, 1]<- rep(1, I*K*R)
  #  X.vec[, 2]<- x1 
  
  # prior for vector beta
  Sigma0<- diag(length(beta.vec))*1000
  Sigma.0.inv = solve(Sigma0)
  beta.0.prior<- c(0,0)
  
  #############################################
  return(list(spatDat = spatDat, y.mat = y.mat, y.vec = y.vec, start.W = start.W, x1=x1,
              # X.mat = X.mat, X.vec=X.vec, 
              X.short = X.short, X.mat = X.mat,
              beta0=beta0, beta1=beta1, Wkeep=Wkeep, 
              p.dat = p.dat,
              sigma2=sigma2, sigma2.s=sigma2.s, Sigma.0.inv = Sigma.0.inv, beta.0.prior=beta.0.prior) )
  
}
