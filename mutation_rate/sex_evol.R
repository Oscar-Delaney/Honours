n = 100 # starting number of individuals of asexual and sexual each
l = 5000 # number of time steps to compute
sd = c(0.01,0.1) # each organisms fitness is initially set to 1 (arbitrary units)
# and upon reproduction, asexual organisms have a standard deviation of 0.01
# and sexuals have a standard deviation of 0.1
cost = 0.02 # sexual reproduction incurs a fitness cost of 0.02

# function to implement reproduction
reproduce <- function(){
  i = sample(seq(1,2*n),1,prob=m[,2]) # choose which organism reproduces
  j = sample(seq(1,2*n)[-i],1) # choose which organism dies (resource constraints)
  # create two new individuals, drawn randomly from a normal distribution
  # depending on the parent, and the mode of reproduction
  new = pmax(rnorm(2,m[i,2]-cost*m[i,1],sd[m[i,1]+1]),0)
  # remove the parent and the deceased organism, and add the two progeny
  m <<- rbind(m[-c(i,j),],matrix(c(rep(m[i,1],2),new),ncol=2))
}

t <- 1:l # vector of timesteps
y <- matrix(0,l,10) # blank matrix to place 10 runs of our simulation in

# create a blank plot to populate as we generate data
plot(NULL,xlim=c(0,l),ylim=c(0,1),xlab="time",ylab="proportion sexual")

for (b in 1:length(y[1,])){
  # initialise a new starting population. Column 1 is the type of organism,
  # 0 means asexual and 1 means sexual, column 2 is that organism's fitness
  m = matrix(c(rep(0, n), rep(1,n*3)), ncol = 2)
  for (a in t) {
    # record the current proportion of sexual organisms
    y[a,b] <- colMeans(m)[1]
    # at each time point, randomly choose an organism to reproduce, weighted by
    # the fitness of each organism
    reproduce()
  }
  lines(y[,b],col=b) # plot that run of the simulation
}

dev.copy(png,"sex_evolution.png", width = 4000, height = 2000, units = "px", res = 300)
dev.off()

