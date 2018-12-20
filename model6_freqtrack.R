# Model 6:  Y - autosome Father's curse - 11 Mar 2018
# Invasion and dynamics of Father's curse Y type with autosomal restorer

# female genotypes  AA   Aa       aa     
#   and fitnesses   1   1-sf/2   1-sf      

#   male genotypes AAY AaY aaY     AAy  Aay   aay
#   and fitnesses   1   1   1      1    1+sm/2  1+sm

# two egg types:     A and  a
# four sperm types:  AY   aY    Ay    ay

library(ggplot2)

generations <- 5000
a <- c(0,1,1,0)
y <- c(0,0,1,1)

set.seed(1865)
plot(a,y,pch=".",main="autosome-Y Father's Curse Dynamics",xlab='a Allele Freq',ylab='y Alelle Freq',font.lab=1.5,cex.main=2,cex.lab=1.35,cex.axis=1.25)
legend(0.8,0.8, legend=c(expression("s" [f]* " < s"[m]), expression("s" [f]* " > s"[m])),
       col=c("red", "black"),pch='.',pt.cex=5)
sel_coef_val.df <- data.frame(matrix(NA, nrow = 100, ncol = 2))
colnames(sel_coef_val.df) <- c('sf_val','sm_val')

for (itest in 1:100){
  #Randomly draw values of sf and sm
  sf_val=.2*runif(1)s
  sm_val=.2*runif(1)
  sel_coef_val.df[itest,] <- c(sf_val,sm_val)
  
  # Specify initial genotype frequencies (3 female, 6 male)
  A = 0.75
  a = 0.25
  Y = 0.75
  y = 0.25
  
  #Initialize genotype frequencies
  #Females
  AA = A * A
  Aa = 2 * A * a
  aa = a * a
  #Males
  AAY = A * A * Y
  AaY = 2 * A * a * Y
  aaY = a * a * Y
  AAy = A * A * y
  Aay = 2 * A * a * y
  aay = a * a * y
  
  for (gen in 1:generations) {
    # Weight each male and female genotype by fitness and normalize
    Fsum=(AA)+(1-(sf_val/2))*(Aa)+(1-sf_val)*(aa)
    AA=AA/Fsum
    Aa=((1-(sf_val/2))*Aa)/Fsum
    aa=((1-sf_val)*aa)/Fsum
    
    Msum=(AAY+AaY+aaY)+(AAy)+(1+(sm_val/2))*(Aay)+(1+sm_val)*(aay)
    AAY=AAY/Msum
    AaY=AaY/Msum
    aaY=aaY/Msum
    AAy=AAy/Msum
    Aay=((1+(sm_val/2))*Aay)/Msum
    aay=((1+sm_val)*aay)/Msum
    
    # Calculate the egg and sperm frequencies
    
    eA=AA+.5*Aa
    ea=aa+.5*Aa
    
    sAY=AAY+.5*AaY
    saY=aaY+.5*AaY
    sAy=AAy+.5*Aay
    say=aay+.5*Aay

    
    # Calculate zygote frequencies the next generation (zygotes are the same freq in both sexes)
    
    AA=eA*(sAY+sAy)
    Aa=eA*(saY+say)+ea*(sAY+sAy)
    aa=ea*(saY+say)
    
    #if (!(isTRUE(all.equal(AA + Aa + aa,1)))) {print ("Fzyg Err")}
    
    AAY=eA*sAY
    AaY=eA*saY+ea*sAY
    aaY=ea*saY
    AAy=eA*sAy
    Aay=eA*say+ea*sAy
    aay=ea*say

    #if (!(isTRUE(all.equal(AAY + AaY + aaY + AAy + Aay+ aay,1)))) {print ("Mzyg Err")}
    
    current_a_female <- (2*aa + Aa)/(2*(AA+Aa+aa))
    current_a_male <- (2*aaY + 2 *aay + AaY + Aay)/(2*(AAY + AaY +aaY + AAy + Aay + aay))
    afreq <- (current_a_female + current_a_male)/2     
    
    yfreq=(AAy+Aay+aay)/(AAY+AaY+aaY+AAy+Aay+aay)
    
    if (sf_val>sm_val){
      points(afreq,yfreq,pch=".")
    } else {
      points(afreq,yfreq,pch=".",col="red")
    } 
  }
}

ggplot(sel_coef_val.df, aes(x=sf_val,y=sm_val)) + geom_point()

