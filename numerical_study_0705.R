########### 이변량정규분포3차원그래프 ##############
if(!require(mvtnorm)) install.packages('mvtnorm'); require(mvtnorm)
if(!require(plotly)) install.packages('plotly'); require(plotly)
if(!require(evd)) install.packages('evd'); require(evd)

set.seed(2019)
n=100
x1 = seq(-5,5,length.out=n)
x2 = seq(-5,5,length.out=n)
mu1 = c(0,-3)
mu2 = c(0,3)
sig = matrix(c(2,1,1,2),nrow=2)
fx = outer(x1, x2, function(x1,x2) { 
  0.4*dmvnorm(cbind(x1,x2),mean=mu1, sigma=sig) +
  0.6*dmvnorm(cbind(x1,x2),mean=mu2, sigma=sig)
  })
plot_ly(x = x1, y = x2, z = fx) %>% add_surface()

# matrix(1:n^2,nrow=n)
ns = 100
s_ind = sample(n^2,ns)
s_row = trunc(s_ind/100)+1 # 1, 100, 101, 10000
s_row[s_ind%%100==0] = s_ind[s_ind%%100==0]/100
s_col = (s_ind/100-s_row+1)*100 
zval = c()
for ( i in 1:ns) { zval[i] =fx[s_row[i],s_col[i]] ;cat(i,fx[s_row[i],s_col[i]],"\n") }
df = data.frame(row=x1[s_row],col=x2[s_col],z=100+zval*100)
plot3d(x=df$row,y=df$col,z=df$z)

xlist = list()
for (i in 1:ns){
xlist[[i]] = rgev(n,loc=df$z[i],scale=40,shape=0.1)
}

