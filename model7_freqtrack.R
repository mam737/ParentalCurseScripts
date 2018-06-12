# Model 7:  Y - X  Father's curse - 12 Mar 2018
# Invasion and dynamics of Father's curse X type with Y-linked restorer

# female genotypes  XX   Xx       xx     
#   and fitnesses   1   1-sf/2   1-sf      

#   male genotypes XY     xY      Xy      xy
#   and fitnesses   1      1      1     1+sm  

# two egg types:     X and  x
# four sperm types:  X   x    Y    y

generations <- 10
x <-c(0,1,1,0)
y <-c(0,0,1,1)
plot(x,y,pch=".",main="X-Y Father's Curse Dynamics",xlab='x Allele Freq',ylab='y Alelle Freq',font.lab=1.5,cex.main=2,cex.lab=1.35,cex.axis=1.25)
legend(0.8,0.8, legend=c(expression("s" [f]* " < s"[m]), expression("s" [f]* " > s"[m])),
       col=c("red", "black"),pch='.',pt.cex=5)

for (itest in 1:100){
  #Randomly draw values of sf and sm
  sf_val=.2*runif(1)
  sm_val=.2*runif(1)
  
  # Specify initial genotype frequencies (3 female, 4 male)
  X = 0.75
  x = 0.25
  Y = 0.75
  y = 0.25
  #Initialize genotype frequencies
  #Females
  XX = X * X
  Xx = 2 * X * x
  xx = x * x
  
  #Males
  XY = X * Y
  xY = x * Y
  Xy = X * y
  xy = x *y
  
  for (gen in 1:generations) {
    # Weight each male and female genotype by fitness and normalize
    Fsum=XX+(1-sf_val/2)*Xx+(1-sf_val)*xx
    XX=XX/Fsum
    Xx=Xx*(1-sf_val/2)/Fsum
    xx=xx*(1-sf_val)/Fsum
    
    Msum=(XY+xY+Xy)+(1+sm_val)*xy
    XY=XY/Msum
    xY=xY/Msum
    Xy=Xy/Msum
    xy=(xy*(1+sm_val))/Msum
    
    # Calculate the egg and sperm frequencies
    
    eX=XX+.5*Xx
    ex=xx+.5*Xx
    if (!(isTRUE(all.equal(eX+ex,1)))) {print ('EggFreqErr')}
    
    sX=.5*(XY + Xy)
    sx=.5*(xY + xy)
    sY=.5*(XY + xY)
    sy=.5*(Xy + xy)
    if (!(isTRUE(all.equal(sX+sx+sY+sy,1)))) {print ('SpermFreqErr')}
    # Calculate zygote frequencies the next generation (zygotes are the same freq in both sexes)
    
    XX=eX*sX
    Xx=eX*sx+ex*sX
    xx=ex*sx
    
    #if (!(isTRUE(all.equal(XX + Xx + xx,1)))) {print (XX + Xx + xx)}
    
    XY=eX*sY
    xY=ex*sY 
    Xy=eX*sy
    xy=ex*sy
    
    #if (!(isTRUE(all.equal(XY + xY + Xy + xy,1)))) {print (XY + xY + Xy + xy)}
    
    current_x_female <- (2*xx + Xx)/(2*(XX+Xx+xx))
    current_x_male <- (xY + xy)/(XY+xY+Xy+xy)
    xfreq <- (current_x_female + current_x_male)/2
    
    yfreq <- (Xy + xy)/(XY + xY + Xy + xy)
    
    if (sf_val>sm_val){
      points(xfreq,yfreq,pch=".")
    } else {
      points(xfreq,yfreq,pch=".",col="red")
    } 
  }
}