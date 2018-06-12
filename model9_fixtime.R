# Model 9:  W-Z nuclear Mother's curse - 13 March 2018
# Invasion and dynamics of nuclear W-linked Mother's curse Z-linked restorer

# female genotypes ZW   Zw    zW    zw  
#   and fitnesses   1   1     1     1+sf  

#   male genotypes  ZZ     Zz      zz
#   and fitnesses   1       1-sm/2       1-sm     

# four egg types:  Z    z     W     w
# two sperm types:  Z  and  z

library(gplots)
library(reshape2)

## 2 paramters of interest, sf and sm
sf_list = seq(0,0.4,length.out = 20)
sm_list = seq(0,0.02,length.out = 20)
generations <- 2500

fix_time.df <- data.frame(matrix(NA,nrow=(length(sf_list)*length(sm_list)), ncol=4))
colnames(fix_time.df) <- c('sf','sm','w_fix_time','z_fix_time')
row_update <- 1

######### MODEL 9 ##########
for (sf_val in sf_list) {
  for (sm_val in sm_list) {
    
    # Specify initial genotype frequencies (4 female, 3 male)
    Z = 0.99
    z = 0.01
    W = 0.99
    w = 0.01
    current_z <- z
    current_w <- w  
    z_fix <- FALSE
    w_fix <- FALSE
    
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
      if (current_w > 0.9999 & w_fix==FALSE) {
        w_fix_gen <- gen
        w_fix <- TRUE
      }
      if (current_z > 0.9999 & z_fix==FALSE) {
        z_fix_gen <- gen
        z_fix <- TRUE
      }
      if ((z_fix == TRUE) & (w_fix == TRUE)) {
        break
      }
      # Weight each male and female genotype by fitness and normalize
      Fsum=ZW+Zw+zW+(zw*(1+sf_val))
      ZW=ZW/Fsum
      Zw=Zw/Fsum
      zW=zW/Fsum
      zw=(zw*(1+sf_val))/Fsum
      
      
      Msum=ZZ+(Zz*(1-(sm_val/2)))+(zz*(1-sm_val))
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
      current_z <- (current_z_female + current_z_male)/2 

      current_w <- (zW + zw)/(ZW + Zw + zW + zw)
    }
    if (w_fix==FALSE) {
      w_fix_gen <- generations
    }
    if (z_fix==FALSE) {
      z_fix_gen <- generations
    }
    fix_time.df[row_update,] <- c(round(sf_val,digits=3),round(sm_val,digits=3),w_fix_gen,z_fix_gen)
    row_update <- row_update + 1
  }
}

pdf('./model9_fix_time_w_heatmap.pdf')
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

pdf('./model9_fix_time_z_heatmap.pdf')
z_fix_time_heatmap <- dcast(fix_time.df[,c(1,2,4)],sm~sf, value.var ='z_fix_time')
row.names(z_fix_time_heatmap) <- z_fix_time_heatmap$sm
z_fix_time_heatmap <- as.matrix(z_fix_time_heatmap[,-1])
if (!all(z_fix_time_heatmap[1,1]==z_fix_time_heatmap)) {
  title = 'Fixation Time for z allele'
  heatmap.2(z_fix_time_heatmap,Rowv=FALSE, Colv=FALSE,dendrogram='none',main=title,xlab='Benefit to Males',ylab='Cost to Females',key=T,key.title='Fixation Time',key.xlab='Generations',density.info='none',trace='none',cexCol=1.25,cexRow=1.25,breaks=seq(0,2500,50),col=redblue(50))
} else {
  print(paste('Fixation Time for z allele is the Same Across all (sf,sm) values: ', z_fix_time_heatmap[1,1]))
}
dev.off()