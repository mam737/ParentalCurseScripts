# Model 8:  W - autosome nuclear Mother's curse - 13 Mar 2018
# Invasion and dynamics of nuclear Mother's curse W type with autosomal restorer

#   female genotypes AAW AaW aaW     AAw  Aaw     aaw
#   and fitnesses    1   1   1       1    1+sf/2  1+sf

#   male genotypes  AA   Aa       aa     
#   and fitnesses   1    1-sm/2   1-sm      

# four egg types:   AW   aW    Aw    aw
# two sperm types:  A  and  a

generations <- 5000
a<-c(0,1,1,0)
w<-c(0,0,1,1)
plot(a,w,pch=".",main="autosome-W Mother's Curse Dynamics",xlab='a Allele Freq',ylab='w Alelle Freq',font.lab=1.5,cex.main=2,cex.lab=1.35,cex.axis=1.25)
legend(0.8,0.9, legend=c(expression("s" [f]* " < s"[m]), expression("s" [f]* " > s"[m])),
       col=c("red", "black"),pch='.',pt.cex=5)

for (itest in 1:100){
  #Randomly draw values of sf and sm
  sf_val=.2*runif(1)
  sm_val=.2*runif(1)
  
  # Specify initial genotype frequencies (6 female, 3 male)
  A = 0.75
  a = 0.25
  W = 0.75
  w = 0.25
  #Initialize genotype frequencies
  #Females
  AAW = A * A * W
  AaW = 2*A * a * W
  aaW = a * a * W
  AAw = A * A * w
  Aaw = 2*A * a * w
  aaw = a * a * w
  
  #Males
  AA = A * A
  Aa = 2*A * a
  aa = a * a
  for (gen in 1:generations) {
    # Weight each male and female genotype by fitness and normalize
    Fsum=AAW+AaW+aaW+AAw+(1+sf_val/2)*Aaw+(1+sf_val)*aaw 
    AAW=AAW/Fsum
    AaW=AaW/Fsum
    aaW=aaW/Fsum
    AAw=(AAw)/Fsum
    Aaw=(Aaw*(1+sf_val/2))/Fsum
    aaw=(aaw*(1+sf_val))/Fsum
    
    Msum=AA+(1-sm_val/2)*(Aa)+(1-sm_val)*(aa)
    AA=AA/Msum
    Aa=Aa*(1-sm_val/2)/Msum
    aa=aa*(1-sm_val)/Msum
    
    # Calculate the egg and sperm frequencies
    
    eAW=AAW+.5*AaW
    eaW=aaW+.5*AaW
    eAw=AAw+.5*Aaw
    eaw=aaw+.5*Aaw
    
    sA=AA+.5*Aa
    sa=aa+.5*Aa
    
    # Calculate zygote frequencies the next generation
    
    AAW=eAW*sA
    AaW=eAW*sa+eaW*sA
    aaW=eaW*sa
    AAw=eAw*sA
    Aaw=eAw*sa+eaw*sA
    aaw=eaw*sa
    
    AA=(eAW + eAw)*sA
    Aa=sA*(eaW+eaw)+sa*(eAW+eAw)
    aa=sa*(eaW+eaw)

    current_a_female <- ((2*aaW) + (2 *aaw) + AaW + Aaw)/(2*(AAW + AaW +aaW + AAw + Aaw + aaw))
    current_a_male <- ((2*aa) + Aa)/(2*(AA+Aa+aa))
    afreq <- (current_a_female + current_a_male)/2      
    
    wfreq <- (AAw + Aaw + aaw)/(AAW+AaW+aaW+AAw+Aaw+aaw)
    
    if (sf_val>sm_val){
      points(afreq,wfreq,pch=".")
    } else {
      points(afreq,wfreq,pch=".",col="red")
    } 
  }
}