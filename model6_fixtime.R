# Model 6:  Y - autosome Father's curse - 11 Mar 2018
# Invasion and dynamics of Father's curse Y type with autosomal restorer

# female genotypes  AA   Aa       aa     
#   and fitnesses   1   1-sf/2   1-sf      

#   male genotypes AAY AaY aaY     AAy  Aay   aay
#   and fitnesses   1   1   1      1    1+sm/2  1+sm

# two egg types:     A and  a
# four sperm types:  AY   aY    Ay    ay


library(gplots)
library(reshape2)

## 2 paramters of interest, sf and sm
sf_list = seq(0,0.02,length.out = 20)
sm_list = seq(0,0.4,length.out = 20)
generations <- 2500

fix_time.df <- data.frame(matrix(NA,nrow=(length(sf_list)*length(sm_list)), ncol=4))
colnames(fix_time.df) <- c('sf','sm','y_fix_time','a_fix_time')
row_update <- 1

######### MODEL 6 ##########
for (sf_val in sf_list) {
  for (sm_val in sm_list) {
    
    # Specify initial genotype frequencies (3 female, 6 male)
    A = 0.99
    a = 0.01
    Y = 0.99
    y = 0.01
    current_a <- a
    current_y <- y  
    a_fix <- FALSE
    y_fix <- FALSE
    
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
      if (current_y > 0.9999 & y_fix==FALSE) {
        y_fix_gen <- gen
        y_fix <- TRUE
      }
      if (current_a > 0.9999 & a_fix==FALSE) {
          a_fix_gen <- gen
          a_fix <- TRUE
      }
      if ((y_fix == TRUE) & (a_fix == TRUE)) {
        break
      }
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
        
      # Calculate zygote frequencies the next generstion (zygotes are the same freq in both sexes)
      
      AA=eA*(sAY+sAy)
      Aa=eA*(saY+say)+ea*(sAY+sAy)
      aa=ea*(saY+say)
      
      #print(paste('F zyg Freq',AA + Aa +aa))
      
      AAY=eA*sAY
      AaY=eA*saY+ea*sAY
      aaY=ea*saY
      AAy=eA*sAy
      Aay=eA*say+ea*sAy
      aay=ea*say
      
      current_a_female <- (2*aa + Aa)/(2*(AA+Aa+aa))
      current_a_male <- (2*aaY + 2 *aay + AaY + Aay)/(2*(AAY + AaY +aaY + AAy + Aay + aay))
      current_a <- (current_a_female + current_a_male)/2     
      
      current_y <- (AAy + Aay + aay)/(AAY+AaY+aaY+AAy+Aay+aay)
      
      }
      if (y_fix==FALSE) {
        y_fix_gen <- generations
      }
      if (a_fix==FALSE) {
        a_fix_gen <- generations
      }
      fix_time.df[row_update,] <- c(round(sf_val,digits=3),round(sm_val,digits=3),y_fix_gen,a_fix_gen)
      row_update <- row_update + 1
  }
}

pdf('./model6_fix_time_y_heatmap.pdf')
y_fix_time_heatmap <- dcast(fix_time.df[,c(1,2,3)],sf~sm, value.var ='y_fix_time')
row.names(y_fix_time_heatmap) <- y_fix_time_heatmap$sf
y_fix_time_heatmap <- as.matrix(y_fix_time_heatmap[,-1])
if (!all(y_fix_time_heatmap[1,1]==y_fix_time_heatmap)) {
  title = 'Fixation Time for y allele'
  heatmap.2(y_fix_time_heatmap,Rowv=FALSE, Colv=FALSE,dendrogram='none',main=title,xlab='Benefit to Males',ylab='Cost to Females',key=T,key.title='Fixation Time',key.xlab='Generations',density.info='none',trace='none',cexCol=1.25,cexRow=1.25,breaks=seq(0,2500,50),col=redblue(50))
} else {
  print(paste('Fixation Time for y allele is the Same Across all (sf,sm) values: ', y_fix_time_heatmap[1,1]))
}
dev.off()

pdf('./model6_fix_time_a_heatmap.pdf')
a_fix_time_heatmap <- dcast(fix_time.df[,c(1,2,4)],sf~sm, value.var ='a_fix_time')
row.names(a_fix_time_heatmap) <- a_fix_time_heatmap$sf
a_fix_time_heatmap <- as.matrix(a_fix_time_heatmap[,-1])
if (!all(a_fix_time_heatmap[1,1]==a_fix_time_heatmap)) {
  title = 'Fixation Time for a allele'
  heatmap.2(a_fix_time_heatmap,Rowv=FALSE, Colv=FALSE,dendrogram='none',main=title,xlab='Benefit to Males',ylab='Cost to Females',key=T,key.title='Fixation Time',key.xlab='Generations',density.info='none',trace='none',cexCol=1.25,cexRow=1.25,breaks=seq(0,2500,50),col=redblue(50))
} else {
  print(paste('Fixation Time for a allele is the Same Across all (sf,sm) values: ', a_fix_time_heatmap[1,1]))
}
dev.off()