# Model 9:  W-Z nuclear Mother's curse - 13 March 2018
# Invasion and dynamics of nuclear W-linked Mother's curse Z-linked restorer

# female genotypes ZW   Zw    zW    zw  
#   and fitnesses   1   1     1     1+sf  

#   male genotypes  ZZ     Zz      zz
#   and fitnesses   1       1-sm/2       1-sm     

# four egg types:  Z    z     W     w
# two sperm types:  Z  and  z

generations <- 10000
z<-c(0,1,1,0)
w<-c(0,0,1,1)
plot(z,w,pch=".",main="W-Z Mother's Curse Dynamics",xlab='z Allele Freq',ylab='w Alelle Freq',font.lab=1.5,cex.main=2,cex.lab=1.35,cex.axis=1.25)
legend(0.8,0.8, legend=c(expression("s" [f]* " < s"[m]), expression("s" [f]* " > s"[m])),
       col=c("red", "black"),pch='.',pt.cex=5)

for (itest in 1:100){
  #Randomly draw values of sf and sm
  sf_val=.2*runif(1)
  sm_val=.2*runif(1)
  
  # Specify initial genotype frequencies (4 female, 3 male)
  Z = 0.75
  z = 0.25
  W = 0.75
  w = 0.25
  
  #Initialize genotype frequencies
  #Females
  ZW = Z * W
  Zw = Z * w 
  zW = z * W
  zw = z * w
  
  #Males
  ZZ = Z * Z
  Zz = 2 * Z * z
  zz = z * z
  
  for (gen in 1:generations) {
    # Weight each male and female genotype by fitness and normalize
    Fsum=ZW+Zw+zW+(zw*(1+sf_val))
    ZW=ZW/Fsum
    Zw=Zw/Fsum
    zW=zW/Fsum
    zw=(zw*(1+sf_val))/Fsum
    
    Msum=ZZ+Zz*(1-sm_val/2)+zz*(1-sm_val)
    ZZ=ZZ/Msum
    Zz=(Zz*(1-sm_val/2))/Msum
    zz=(zz*(1-sm_val))/Msum
    
    # Calculate the egg and sperm frequencies
    
    eZ=.5*(ZW+Zw)
    ez=.5*(zW+zw)
    eW=.5*(ZW+zW)
    ew=.5*(Zw+zw)
    
    sZ=ZZ+.5*Zz
    sz=zz+.5*Zz
    
    # Calculate zygote frequencies the next generstion
    
    ZW=eW*sZ
    Zw=ew*sZ 
    zW=eW*sz
    zw=ew*sz

    
    ZZ=eZ*sZ
    Zz=eZ*sz+ez*sZ
    zz=ez*sz
    
    current_z_female <- (zW + zw)/(ZW + Zw + zW + zw) 
    current_z_male <- (Zz + 2*zz)/(2*(ZZ + Zz + zz))
    zfreq <- (current_z_female + current_z_male)/2      
    
    wfreq <- (Zw + zw)/(ZW + Zw + zW + zw)

    
    if (sf_val>sm_val){
      points(zfreq,wfreq,pch=".")
    } else {
      points(zfreq,wfreq,pch=".",col="red")
    } 
    
  }
}