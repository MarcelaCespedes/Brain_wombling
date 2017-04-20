########################################################################
# started: Thursday 3nd Nov
# 35 region MH within Gibbs for unbalanced spatial data set
# with option to initialise chains at solution or away from solution
#

#-------------------------------------------------------

library(ggplot2)
library(MASS)
library(stats) # with 1/gamma(1, shape, rate)
library(raster)
library(mvtnorm)
library(mefa)

rm(list = ls())

source("IRK_35Unbal_genData.r")
pdf.name<- "35Unbal_high_rhos2s.pdf"

K=35
I=30 
R.vec<- c(4, 3, 5, 2, 4, 3, 4, 5, 4, 3,
          4, 3, 5, 2, 4, 3, 4, 5, 4, 3, 
          4, 3, 5, 2, 4, 3, 4, 5, 4, 3)  # <=== this vector to be of same length as I
length(R.vec)

##
## set the number of runs for MCMC scheme
##

mcmc = 10000

no.female=round(0.5*I); 

##
## Set values to generate data
##
rho = 0.9

no.Obs<- sum(R.vec)*K  # <=== no obs
no.Obs

s2s =  2   # values set as per the manuscripts simulation study    
s2= 0.5       
b0 = 3
b1 = 0.5

# Generate data from model
dat<- IRK_35Unbal_genData(I = I, no.female=no.female, s2s=s2s, s2=s2, b0=b0, b1=b1, rho.val= rho, R.vec= R.vec,
                        mcmc = mcmc)

# view data
x11()
plot(dat$p.dat)

####################################################################################
# Set up chains for wobmling algorithm

beta0= dat$beta0
beta0[1] = runif(1, min = 2, max = 3) #rnorm(1) 

beta1= dat$beta1      # <=== initialised empty chains for model
beta1[1] = runif(1, min = -0.5, max = 0.5)  #rnorm(1)

Wkeep = dat$Wkeep   # <== To start at the solution

load("contig_35Reg.Rdata")
Wkeep[,,1] = contig_35Reg  # <=== initialise off solution

cNt<- matrix(0, K, K)

sigma2 = dat$sigma2
sigma2[1] = runif(1, min = 0.0001, max = 1)

sigma2.s = dat$sigma2.s
sigma2.s[1] =  runif(1, min = 0.0001, max = 3)

y.vec = dat$y.vec  # <===== data 
y.mat = dat$y.mat
X.mat = dat$X.mat

X.mat.long = as.matrix(rep(data.frame(X.mat), (K*R.vec)))

Sigma.0.inv = dat$Sigma.0.inv   # <======== beta priors
mean.0 = c(0, 0)

# I think I'll have to also specify the prior for s2, s2s, etc
# and finadle with X.mat, etc

# prior for s2s
sh.s2s = 1; rat.s2s = 2
#dinvgamma(s2s, shape = sh.s2s, rate = rat.s2s)  # check support for prior, want to have high density support

# prior for s2
sh.s2 = 1; rat.s2 = 0.5
#dinvgamma(s2, shape = sh.s2, rate = rat.s2) # check support for prior, want to have high density support

#####
# priors for beta.vec
Sigma.o.inv<- matrix(c(1,0,0,1), 2, 2) + diag(2)*0.5
Sigma.o.inv
mu.o<- c(beta0[1], beta1[1]) #rep(1, 2)

#####
# set random effects

b.mat<- matrix(0, nrow = I, ncol = K)
# sum over all the replicates for person i - note this doesn't change... can do this external to the code
sum.i.rep<- matrix(0, nrow = I, ncol = K)

for(i in 1:I){
  
  if(i == 1){
    sum.i.rep[i,] <- colSums(y.mat[ (i + (i-1)*R.vec[i]) :R.vec[i], ]) # checked :)
  }else{
    sum.i.rep[i, ]<- colSums(y.mat[ (sum(R.vec[1:i-1])+1)  :sum(R.vec[1:i]), ])
  }
}

e.1<- rep(1, K)

########################################################################################
# start Gibbs
########################################################################################

#----------------------------------------------------------------------
tiMe<- proc.time() 
pb <- txtProgressBar(style = 3)

for(t in 2:mcmc){
  
  ##################################################################################
  # draw the random effects - this will be used for the other samplers below
  # b_i
  
  beta.vec <- c(beta0[t-1], beta1[t-1])
  Q.inv = rho*(diag(K)*rowSums(Wkeep[,,t-1]) - Wkeep[,,t-1]) + (1 - rho)*diag(K) 
  
  for(i in 1:I){
    Omega = solve(R.vec[i]/sigma2[t-1]*diag(K) + 1/sigma2.s[t-1]*Q.inv)
    
    mu = Omega%*%(1/sigma2[t-1]*sum.i.rep[i,] - R.vec[i]/sigma2[t-1]*X.mat[i,]%*%beta.vec*e.1)
    b.mat[i, ]<- rmvnorm(1, mean = mu, sigma = Omega)
  }  
  
  #round(b.mat, 2)  # compare random effects with solution
  #round(dat$spatDat, 2)
  
  ###############################################
  # update sigma2.s
  #sigma2.s[t] <- sigma2.s[t-1]
  
  shape.s2s = (I*K + 2*sh.s2s)/2    
  #shape.s2s                         
  
  # need to add over all the random effects
  sum.rand.eff<- c()
  for(ii in 1:I){
    sum.rand.eff[ii]<- b.mat[ii,]%*%Q.inv%*%b.mat[ii,]
  }
  rate.s2s = 0.5*(sum(sum.rand.eff)) + rat.s2s
  rate.s2s
  
  sigma2.s[t]<- 1/rgamma(1, shape = shape.s2s, rate = rate.s2s)
  head(sigma2.s)

  
  ################################################
  # update W
  #Wkeep[,,t] <- Wkeep[,,t-1]
  
  ck.prob<- ck.proposal.m<- matrix(0, K, K)  # Try to do a Metropolis Hastings
  w.curr<- Wkeep[,,t-1]
  for(row in 1:(K-1)){
    for(col in (row+1):K){
      
      w.star<- w.curr
      w.star[row, col] <- w.star[col, row]<- ck.proposal.m[row, col]<- rbinom(1, 1, prob = 0.5)  # prior for each w_kj
      #w.star
      
      Q.curr = rho*(diag(K)*rowSums(w.curr) - w.curr) + (1 - rho)*diag(K) 
      #        Q.curr
      Q.star =rho*(diag(K)*rowSums(w.star) - w.star) + (1 - rho)*diag(K) 
      #        Q.star
      
      #### check notes and add the determinant of Q.star and Q.curr
      prop.star<- curr<- rep(0, I)
      
      for(i in 1:I){
        prop.star[i] = 0.5*log(det(Q.star)) -1/(2*sigma2.s[t-1])*t(b.mat[i,])%*%Q.star%*%b.mat[i,]
        curr[i] = 0.5*log(det(Q.curr))-1/(2*sigma2.s[t-1])*t(b.mat[i,])%*%Q.curr%*%b.mat[i,]
      }
      
      acc.rej = exp(sum(prop.star) - sum(curr) ) 
      #acc.rej
      ck.prob[row, col]<- acc.rej
      
      if(acc.rej > runif(1)){
        w.curr = w.star
        cNt[row, col]<- cNt[col, row]<- cNt[row, col] + 1}
    }
  }
  #w.curr
  Wkeep[,,t]<-w.curr
  
  ################################################
  # update sigma2
  #sigma2[t] <- sigma2[t-1]
  
  shape.s2 = no.Obs/2 + sh.s2  # note the shape for this distribution won't change
  shape.s2
  # convert everything to vector form
  rand.long<- as.matrix(rep(data.frame(b.mat), R.vec))# repeat every random effects according to R.vec
  rand.long2<- as.vector(t(rand.long))
  
  sum.lin.pred<- c()
  for(all in 1:(no.Obs)){
    sum.lin.pred[all]<- (y.vec[all] - X.mat.long[all,]%*%beta.vec - rand.long2[all])^2
  }
  rate.s2 = sum(sum.lin.pred)/2 + rat.s2
  rate.s2
  sigma2[t]<- 1/rgamma(1, shape = shape.s2, rate = rate.s2)
  
  #head(sigma2)
  
  # note: From the form of the shape and scale GIbbs samplers, we can see that it is robust to the
  # choice of prior: The prior value has very little effect on the shape and rate of the posterior of s2 
  ################################################
  # update beta's
  #beta0[t] <- beta0[t-1]
  #beta1[t] <- beta1[t-1]
  
  Omega.b<- solve(t(X.mat.long)%*%X.mat.long/sigma2[t] + Sigma.o.inv)
  
  #b.rep<- as.matrix(rep(data.frame(b.mat), each = R)) # repeat each b.mat column R times
  #b.vec<- as.vector(t(b.rep))
  
  mu = t(t(X.mat.long)%*%(y.vec - rand.long2)/sigma2[t] +  t(mu.o%*%Sigma.o.inv))%*%Omega.b
  beta.vec<- rmvnorm(1, mean=mu, sigma=Omega.b)
  
  beta0[t]<- beta.vec[1]
  beta1[t]<- beta.vec[2]
  
  
  ####################################################################
  prog<- t/mcmc 
  setTxtProgressBar(pb, prog) 
}

(proc.time() - tiMe)/60 # time in minutes


#####################   density -----------------------------------------------
p.s2s<- ggplot(data.frame(x = sigma2.s), aes(x=x)) + geom_density() + geom_vline(xintercept = s2s, colour = "red") + 
  theme_bw() + theme(legend.position = "none") + ggtitle("s2s")

p.s2<- ggplot(data.frame(x = sigma2), aes(x=x)) + geom_density() + geom_vline(xintercept = s2, colour = "red") + 
  theme_bw() + theme(legend.position = "none") + ggtitle("s2")

p.b0<- ggplot(data.frame(x = beta0), aes(x=x)) + geom_density() + geom_vline(xintercept = b0, colour = "red") + 
  theme_bw() + theme(legend.position = "none") + ggtitle("b0")

p.b1<- ggplot(data.frame(x = beta1), aes(x=x)) + geom_density() + geom_vline(xintercept = b1, colour = "red") + 
  theme_bw() + theme(legend.position = "none") + ggtitle("b1")

source("multiplot.r")
x11()
multiplot(p.s2s, p.s2, p.b0, p.b1, cols = 2)

#################### trace -----------------------------------------------------
x = seq(1:mcmc)
trace.s2s<- ggplot(data.frame(x=x, y = sigma2.s), aes(x=x, y=y)) + geom_line() + geom_hline(yintercept = s2s, colour = "red") + 
  theme_bw() + theme(legend.position = "none") + ggtitle("s2s")

trace.s2<- ggplot(data.frame(x=x, y = sigma2), aes(x=x, y=y)) + geom_line() + geom_hline(yintercept = s2, colour = "red") + 
  theme_bw() + theme(legend.position = "none") + ggtitle("s2")

trace.b0<- ggplot(data.frame(x = x, y = beta0), aes(x=x, y=y)) + geom_line() + geom_hline(yintercept = b0, colour = "red") + 
  theme_bw() + theme(legend.position = "none") + ggtitle("b0")

trace.b1<- ggplot(data.frame(x = x, y=beta1), aes(x=x,y=y)) + geom_line() + geom_hline(yintercept = b1, colour = "red") + 
  theme_bw() + theme(legend.position = "none") + ggtitle("b1")

#source("multiplot.r")
x11()
multiplot(trace.s2s, trace.s2, trace.b0, trace.b1, cols = 2)

################# plot W
no.off.diag<- K*(K-1)/2
W.param<- matrix(0, nrow = dim(Wkeep)[3], ncol = no.off.diag) # take the off-diagonals of the matrices pkeep and Wkeep
couNt = 1

for(row in 1:(K-1) ){  # loop thru upper triangular
  for(col in (row + 1): K){
    W.param[, couNt]<- Wkeep[row, col,]
    couNt = couNt + 1
  }
}

w.m<- matrix(0, K, K); count<- 1
w.post.med<- apply(W.param, 2, FUN = mean)

for(row in 1:(K-1) ){ 
  for(col in (row + 1): K){
    w.m[row, col]<- w.post.med[count]
    w.m[col,row]<- w.m[row, col]
    count<- count + 1
  }
}

r.w <- r.acc.rate<- r.start<-raster(xmn = 0, xmx =K, ymn = 0, ymx = K, nrows = K, ncols = K)
r.w[]<- as.vector(w.m) 
r.acc.rate[]<- as.vector(cNt/mcmc)
r.start[]<- as.vector(Wkeep[,,1])

x11()
par(mfrow = c(1,3))
plot(r.w, main = "W posterior mean")
plot(r.acc.rate, main = "Acceptance rate for W")
plot(r.start, main = "Start W",legend = F)


######################################################
# now apply thinning and burn in
#source("MCMC_35_diags.r")

burnIn = floor(mcmc/5)
thin = 5

op<- MCMC_35_diags(Wkeep = Wkeep, thin=thin, burnIn=burnIn, beta0=beta0, beta1=beta1,
                   b0.sol = b0, b1.sol =b1, s2=s2, s2s=s2s,
                   sigma2=sigma2, sigma2.s=sigma2.s, mcmc=mcmc, accept.prob = cNt/mcmc)

op$length.chain # length of remaining chain

source("multiplot.r")
op$b0.diag
windows()
multiplot(op$b0.ac, op$b0.trace, op$b0.dens, cols =3)

op$b1.diag
windows()
multiplot(op$b1.ac, op$b1.trace, op$b1.dens, cols =3)

op$sigma2.diag
windows()
multiplot(op$sigma.ac, op$sigma.trace, op$sigma.dens, cols =3)

op$ss2.diag
windows()
multiplot(op$ss.ac, op$ss.trace, op$ss.dens, cols =3)

windows()
par(mfrow = c(1,2))
plot(op$r.w, main = "Posterior mean W")
plot(op$r.p, main = "Acceptance ratio")

###################################################################################
###################################################################################
# save to pdf
pdf.name

pdf(pdf.name, onefile = T)
#windows()
#plot(1:20, type = 'n', bty = "n", xlab = "", ylab = "")
plot(1:20, type = 'n', xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
text(10, 16, paste("rho = 0.99"))
text(10, 14, paste("I ppl: ", I, sep = ""))
text(10, 13, paste("R unbalanced:\n ", R.vec, sep = ""))
text(10, 12, paste("K regions: ", K, sep = ""))
text(10, 11, paste("No mcmc: ", mcmc, sep = ""))

text(10, 10, paste("sol.b0:", b0, 
                   "     low.95.ci:", round(quantile(beta0, probs = 0.025),2), 
                   "     post.mean:", round(mean(beta0), 2), 
                   "     high.95.ci:", round(quantile(beta0, probs = 0.975),2), sep = "" ))

text(10, 8, paste("sol.b1:", b1, 
                  "     low.95.ci:", round(quantile(beta1, probs = 0.025),2), 
                  "     post.mean:",round(mean(beta1), 2),  
                  "     high.95.ci:",round(quantile(beta1, probs = 0.975),2), sep = ""))

text(10, 6, paste("sol.s2:", s2, 
                  "     low.95.ci:", round(quantile(sigma2, probs = 0.025),2), 
                  "     post.mean:", round(mean(sigma2), 2), 
                  "     high.95.ci:", round(quantile(sigma2, probs = 0.975),2), sep = ""))

text(10, 4, paste("sol.s2s:", s2s, 
                  "     low.95.ci:", round(quantile(sigma2.s, probs = 0.025),2), 
                  "     post.mean:", round(mean(sigma2.s), 2),  
                  "     high.95.ci:",round(quantile(sigma2.s, probs = 0.975),2), sep = ""))
plot(dat$p.dat)
multiplot(p.s2s, p.s2, p.b0, p.b1, cols = 2)
multiplot(trace.s2s, trace.s2, trace.b0, trace.b1, cols = 2)
plot(r.w, main = "W posterior mean")
plot(r.acc.rate, main = "Acceptance rate for W")
plot(r.start, main = "Starting W", legend = F)
multiplot(op$b0.ac, op$b0.trace, op$b0.dens, cols =3)
multiplot(op$b1.ac, op$b1.trace, op$b1.dens, cols =3)
multiplot(op$sigma.ac, op$sigma.trace, op$sigma.dens, cols =3)
multiplot(op$ss.ac, op$ss.trace, op$ss.dens, cols =3)
dev.off()

