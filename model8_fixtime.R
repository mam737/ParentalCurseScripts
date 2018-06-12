# Model 8:  W - autosome nuclear Mother's curse - 13 Mar 2018
# Invasion and dynamics of nuclear Mother's curse W type with autosomal restorer

#   female genotypes AAW AaW aaW     AAw  Aaw     aaw
#   and fitnesses    1   1   1       1    1+sf/2  1+sf

#   male genotypes  AA   Aa       aa     
#   and fitnesses   1    1-sm/2   1-sm      

# four egg types:   AW   aW    Aw    aw
# two sperm types:  A  and  a

library(gplots)
library(reshape2)

## 2 paramters of interest, sf and sm
sf_list = seq(0,0.4,length.out = 20)
sm_list = seq(0,0.02,length.out = 20)
generations <- 2500

fix_time.df <- data.frame(matrix(NA,nrow=(length(sf_list)*length(sm_list)), ncol=4))
colnames(fix_time.df) <- c('sf','sm','w_fix_time','a_fix_time')
row_update <- 1

######### MODEL 8 ##########
for (sf_val in sf_list) {
  for (sm_val in sm_list) {
    
    # Specify initial genotype frequencies (6 female, 3 male)
    A = 0.99
    a = 0.01
    W = 0.99
    w = 0.01
    current_a <- a
    current_w <- w  
    a_fix <- FALSE
    w_fix <- FALSE
    
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
      if (current_w > 0.9999 & w_fix==FALSE) {
        w_fix_gen <- gen
        w_fix <- TRUE
      }
      if (current_a > 0.9999 & a_fix==FALSE) {
        a_fix_gen <- gen
        a_fix <- TRUE
      }
      if ((w_fix == TRUE) & (a_fix == TRUE)) {
        break
      }
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
      
      # Calculate zygote frequencies the next generstion
      
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
      current_a <- (current_a_female + current_a_male)/2   
      
      current_w <- (AAw + Aaw + aaw)/(AAW+AaW+aaW+AAw+Aaw+aaw)
    }
    if (w_fix==FALSE) {
      w_fix_gen <- generations
    }
    if (a_fix==FALSE) {
      a_fix_gen <- generations
    }
    fix_time.df[row_update,] <- c(round(sf_val,digits=3),round(sm_val,digits=3),w_fix_gen,a_fix_gen)
    row_update <- row_update + 1
  }
}

pdf('./model8_fix_time_w_heatmap.pdf')
w_fix_time_heatmap <- dcast(fix_time.df[,c(1,2,3)],sm~sf, value.var ='w_fix_time')
row.names(w_fix_time_heatmap) <- w_fix_time_heatmap$sm
w_fix_time_heatmap <- as.matrix(w_fix_time_heatmap[,-1])
if (!all(w_fix_time_heatmap[1,1]==w_fix_time_heatmap)) {
  title = 'Fixation Time for w allele'
  heatmap.2(w_fix_time_heatmap,Rowv=FALSE, Colv=FALSE,dendrogram='none',main=title,xlab='Benefit to Males',ylab='Cost to Females',key=T,key.title='Fixation Time',key.xlab='Generations',density.info='none',trace='none',cexCol=1.25,cexRow=1.25,breaks=seq(0,2500,50),col=redblue(50))
} else {
  print(paste('Fixation Time for w allele is the Same Across all (sf,sm) values: ', w_fix_time_heatmap[1,1]))
}
dev.off()

pdf('./model8_fix_time_a_heatmap.pdf')
a_fix_time_heatmap <- dcast(fix_time.df[,c(1,2,4)],sm~sf, value.var ='a_fix_time')
row.names(a_fix_time_heatmap) <- a_fix_time_heatmap$sm
a_fix_time_heatmap <- as.matrix(a_fix_time_heatmap[,-1])
if (!all(a_fix_time_heatmap[1,1]==a_fix_time_heatmap)) {
  title = 'Fixation Time for a allele'
  heatmap.2(a_fix_time_heatmap,Rowv=FALSE, Colv=FALSE,dendrogram='none',main=title,xlab='Benefit to Males',ylab='Cost to Females',key=T,key.title='Fixation Time',key.xlab='Generations',density.info='none',trace='none',cexCol=1.25,cexRow=1.25,breaks=seq(0,2500,50),col=redblue(50))
} else {
  print(paste('Fixation Time for a allele is the Same Across all (sf,sm) values: ', a_fix_time_heatmap[1,1]))
}
dev.off()