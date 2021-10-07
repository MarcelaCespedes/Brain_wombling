###################################################
# Generate data from model for Wombling algorithm
# MRI cortical thickness data generated for 35 regions of interest (ROI) of the left hemisphere of the
# brain. ROIs were parcellated according to the Automated Anatomical Labeling (AAL) atlas. See book chapter for listing of ROIs (https://www.springer.com/gp/book/9783030425524)

# Repeated measure MRI ROI data generated in an unbalanced design. This means simulated participants did not have the same number of repeated measures (replicates) to better represent real-life scenario. 

# Spatial association between the ROIs is simulated by a Conditional Autoregressive (CAR) spatial random effect. 

# Covariates measures across participants only inlcude Sex effects in this simulated data set.

# I: number of participants
# K: total numberof ROIs
# R.vec: vector with the number of repeated MRI ROI measures per person, length(R.vec) = I



# Marcela Cespedes
# 7 October 2021

Unbalanced_genData<- function(I = 5, # number of patients 
                         no.female = 3, # number of females, Males = I-no.female 
                         s2s = 1, # spatial scale variance term 
                         R.vec = rep(2, 5), # vector of repeated measures across each person 
                         s2 = 0.5, # residual variance 
                         b0 = 3, # average cortical thickness across all ROIs 
                         b1=0.5, # Sex effect 
                         rho.val=0.9){ # level of spatial correlation 
                         
  # default set total number of ROI to 35
  K<- 35
  no.obs<- sum(R.vec)*K
  #no.obs # check
  
  #######################################################################
  # simulate participant sex values
  # binary variable, female = 1, male is baseline = 0
  x1 = c(rep(1, K*sum(R.vec[1:no.female])), 
         rep(0,K*sum(R.vec[(1 + no.female):I]) )   )
  x1.short = c(rep(1, no.female), rep(0, I-no.female))
  
  
  ########################################################################
  # simulate spatial random effects
  
  # load("contig_35Reg.Rdata") # use to generate data according to contiguous matrix
  load("strong_spat35Reg.Rdata") # leading diagonal correlated ROIs
  W = strong_spat35Reg
  diag(W)<- 0  
  #W # check
  
  start.W<- W
  Q = rho.val*(diag(K)*rowSums(W) - W) + (1-rho.val)*diag(K)
  # generate spatial random effects across I participants
  spatDat<- mvrnorm(I, mu = rep(0, K), Sigma = s2s*solve(Q) ) #<- spatial random effects simulated 
  
  library(mefa)
  vec.d<- as.matrix(rep(data.frame(spatDat), R.vec))# repeat every random effect according to R.vec elements
  vec.data<- as.vector(t(vec.d))
  
  
  ######################################################################
  # simulate ROI MRI cortical thickness data (outcome), taking into account Sex effects
  y = b0 + b1*x1 + vec.data + rnorm(no.obs, mean = 0, sd = sqrt(s2))
  y.vec = y # <-- outcome simulated
  
  no.obs.ppl<- R.vec*K
  dat<- data.frame(y= y.vec, region = factor(rep(1:K, times = sum(R.vec)*K  )), person = factor(rep(1:I, no.obs.ppl) ) )
  
  # check simulated cortical thickness measures across K ROIs for I participants
  #x11()
  p.dat<- ggplot(dat, aes(x=region, y=y, colour = region)) + geom_point()+
    ggtitle(paste("35 region CAR ROI MRI simulated data for ", I, "participants with unlabanced replicates", sep = "")) + theme_bw() + facet_wrap(~person, ncol = 5)
  #p.dat
  
  
  
  
  
  #####################################################################
  #####################################################################
  ## Set up data frames for MCMC sampling
  ##
  
  # convert longitudinal outcome to wide format. Each column corresponds to ROI
  tR<- matrix(y, nrow = K, ncol = sum(R.vec))
  y.mat<- t(tR)
  
  # W initialised at solution
  Wkeep<-array(0, dim = c(K,K, mcmc)); Wkeep[,,1]<- start.W;   
  
  # fixed effects and variance terms initialised at solution
  beta0 <- beta1<- sigma2<- sigma2.s<- rep(0, mcmc) 
  
  beta0[1] = b0; beta1[1] = b1; 
  sigma2[1] = s2; sigma2.s[1] = s2s; beta.vec = c(b0,b1)
  
  # need to have data in both matrix and vector form
  X.mat = matrix(0, nrow = I, ncol = 2)
  X.mat[, 1]<- rep(1, I)
  X.mat[, 2]<- x1.short
  
  #############################################
  #############################################
  return(list(#spatDat = spatDat, 
              y.mat = y.mat, 
              y.vec = y.vec, 
              start.W = start.W, x1=x1,
              X.mat = X.mat,
              beta0=beta0, beta1=beta1, Wkeep=Wkeep, 
              p.dat = p.dat,
              sigma2=sigma2, sigma2.s=sigma2.s) )
  
}
