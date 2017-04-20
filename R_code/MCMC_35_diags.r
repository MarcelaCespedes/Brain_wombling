###############################################################
# process the Gibbs sampler 

MCMC_35_diags<- function(Wkeep = Wkeep, accept.prob =accept.prob,
                          beta0 = beta0, b0.sol = b0.sol, beta1 = beta1, b1.sol = b1.sol,s2 = s2, s2s =s2s,
                          burnIn = 3000, thin=2, mcmc =10000, sigma2= sigma2, sigma2.s = sigma2.s){
  
  library(ggplot2)
  library(reshape2)
  library(raster)
  
  w.m<- Wkeep[,,burnIn:mcmc]
  
  N<- 35; couNt = 1
  no.off.diag<- N*(N-1)/2
  
  W.param<- matrix(0, nrow = dim(w.m)[3], ncol = no.off.diag) # take the off-diagonals of the matrices pkeep and Wkeep
  for(row in 1:(N-1) ){  # loop thru upper triangular
    for(col in (row + 1): N){
      W.param[, couNt]<- w.m[row, col,]
      couNt = couNt + 1
    }
  }
  
  n.delete<- function(dataframe, n){ # n=1 == no thinning
    out<- dataframe[seq(n, to=nrow(dataframe), by=n),]
    return(out)
  }
  
  thin.val<-thin
  
  w.mat<- as.data.frame(W.param)
  w.p<- n.delete(dataframe = w.mat, n=thin.val)
  
  #####################################################################################
  # plot neighbourhood matrix
  w.mAT<-matrix(0, N, N); count<- 1
  w.post.med<- apply(w.p, 2, FUN = mean) # raster is computed for non-thinned W (?)
  
  for(row in 1:(N-1) ){ 
    for(col in (row + 1): N){
      w.mAT[row, col]<- w.post.med[count]
      w.mAT[col,row]<- w.mAT[row, col]
      count<- count + 1
    }
  }
  
  r.w <- r.p<- raster(xmn = 0, xmx =N, ymn = 0, ymx = N, nrows = N, ncols = N)
  r.w[]<- as.vector(w.mAT) 
  
  r.p[]<- as.vector(accept.prob)
  
  # windows()
  #  plot(r.w, legend = TRUE, main ="Posterior median W"  )
  #####################################################################################
  # beta0 diagnostics  
  
  b0<- beta0[burnIn:mcmc]
  b0<- n.delete(dataframe=data.frame(b0), n=thin.val) 
  
  b0.diag = data.frame(actual.B0 = b0.sol, low.95.ci = quantile(b0, probs = 0.025),post.mean = mean(b0),  
                       high.95.ci = quantile(b0, probs = 0.975))
  
  # AUTOCORRELATION, for each parameter estimate
  #windows()
  b0.auto.corr<- with(acf(b0, plot=FALSE), data.frame(lag, acf))
  
  b0.ac <- ggplot(data = b0.auto.corr, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) + theme_bw() +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Autocorr for b0")
  
  
  th3<- data.frame(b0 = b0, id = seq(1:length(b0)))
  #windows()
  b0.trace<- ggplot(data = th3, aes(x=id, y=b0)) + geom_line() +
    theme(legend.position="top") +  labs(x="MCMC iteration", y="Simulation from b0 marginal posterior") +
    ggtitle("Trace plot for b0") + theme_bw()
  
  b0.dens<- ggplot(data = th3, aes(x=b0)) + geom_density() + 
    ggtitle("Density for b0 marginal ") + theme_bw()
  
  ###################################################################################################
  # beta1 diagnostics  
  
  b1<- beta1[burnIn:mcmc]
  b1<- n.delete(dataframe=data.frame(b1), n=thin.val) 
  
  b1.diag = data.frame(actual.B1 = b1.sol, low.95.ci = quantile(b1, probs = 0.025),post.mean = mean(b1),  
                       high.95.ci = quantile(b1, probs = 0.975))
  
  # AUTOCORRELATION, for each parameter estimate
  #windows()
  b1.auto.corr<- with(acf(b1, plot=FALSE), data.frame(lag, acf))
  
  b1.ac <- ggplot(data = b1.auto.corr, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) + theme_bw() +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Autocorr for b1")
  
  
  th3<- data.frame(b1 = b1, id = seq(1:length(b1)))
  #windows()
  b1.trace<- ggplot(data = th3, aes(x=id, y=b1)) + geom_line() +
    theme(legend.position="top") +  labs(x="MCMC iteration", y="Simulation from b1 marginal posterior") +
    ggtitle("Trace plot for b1") + theme_bw()
  
  b1.dens<- ggplot(data = th3, aes(x=b1)) + geom_density() + 
    ggtitle("Density for b1 marginal ") + theme_bw()
  
  ######################################################################################
  # residual variance
  #head(sigmakeep)
  s<- sigma2[burnIn:mcmc]
  s<- n.delete(dataframe=data.frame(s), n=thin.val) 
    
  sigma2.diag<- data.frame(s2.sol = s2, low.95.ci = quantile(s, probs = 0.025),post.median = median(s),  
                           high.95.ci = quantile(s, probs = 0.975))
    
  sigma.auto.corr<- with(acf(s, plot=FALSE), data.frame(lag, acf))
    
  sigma.ac <- ggplot(data = sigma.auto.corr, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) + theme_bw() +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Autocorr for s")
    
  th3<- data.frame(s = s, id = seq(1:length(s)))
    #windows()
  sigma.trace<- ggplot(data = th3, aes(x=id, y=s)) + geom_line() +
    theme(legend.position="top") +  labs(x="MCMC iteration", y="Simulation from s2 marginal posterior") +
    ggtitle("Trace plot for s2") + theme_bw()
    
  sigma.dens<- ggplot(data = th3, aes(x=s)) + geom_density() + 
    ggtitle("Density for s2 marginal ") + theme_bw()
  
  ######################################################################################
  # note: we also estimate sigma.s <= within group variance
  #head(sigma.skeep)
  
  ss<- sigma2.s[burnIn:mcmc]
  ss<- n.delete(dataframe=data.frame(ss), n=thin.val) 
    
    #length(ss)
  ss2.diag= data.frame(s2s.sol = s2s, low.95.ci = quantile(ss, probs = 0.025),post.median = median(ss),  
                       high.95.ci = quantile(ss, probs = 0.975))
    
  ss.auto.corr<- with(acf(ss, plot=FALSE), data.frame(lag, acf))
    
  ss.ac <- ggplot(data = ss.auto.corr, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) + theme_bw() +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + ggtitle("Autocorr for s2s")
    
  th3<- data.frame(ss = ss, id = seq(1:length(ss)))
    #windows()
  ss.trace<- ggplot(data = th3, aes(x=id, y=ss)) + geom_line() +
    theme(legend.position="top") +  labs(x="MCMC iteration", y="Simulation from s2s marginal posterior") +
    ggtitle("Trace plot for s2s") + theme_bw()
    
  ss.dens<- ggplot(data = th3, aes(x=ss)) + geom_density() + 
    ggtitle("Density for ss marginal ") + theme_bw()
    
    ############################################################################################################
  return(list(length.chain = dim(w.p)[1], 
              W.chains = w.p,b0=b0, b1=b1, w.mAT = w.mAT,
              r.w= r.w,  r.p= r.p,
              b0.diag = b0.diag, b1.diag = b1.diag, 
              b0.ac = b0.ac, b0.trace = b0.trace, b0.dens = b0.dens,
              b1.ac = b1.ac, b1.trace = b1.trace, b1.dens = b1.dens,
              b0.chains = b0, b1.chains = b1, sigma2.chains = s, ss.chains = ss,
              sigma2.diag= sigma2.diag, sigma.ac = sigma.ac, sigma.trace = sigma.trace, sigma.dens=sigma.dens,
              ss2.diag = ss2.diag, ss.ac = ss.ac, ss.trace=ss.trace, ss.dens=ss.dens))
} # end function
