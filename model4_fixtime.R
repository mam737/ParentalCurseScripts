# Model 4:  mito-Z Mother's curse - 11 March 2018
# Invasion and dynamics of Mother's curse mitochondrial type with Z-linked restorer

# female genotypes ZWM   zWM    ZWm       zWm
#   and fitnesses   1    1      1+sf      1+sf  

#   male genotypes  ZZM  ZzM    zzM   ZZm     Zzm             zzm
#   and fitnesses    1    1      1     1-sm    1-sm +sz/2     1-sm + sz

# six egg types:   ZM   Zm   zM   zm    WM   Wm    
# two sperm types:  Z  and   z

# Specify initial genotype frequencies (4 female, 6 male)


# Weight each male and female genotype by fitness and normalize

library(gplots)
library(reshape2)

## 3 paramters of interest, sf, sm, and sa
sf_list = seq(0,0.4,length.out = 20)
sm_list = seq(0,0.2,length.out = 20)
sz_list = seq(0,0.2,length.out = 20)
generations <- 1000

fix_time.df <- data.frame(matrix(NA,nrow=(length(sf_list)*length(sm_list)*length(sz_list)), ncol=5))
colnames(fix_time.df) <- c('sf','sm','sz' ,'z_fix_time','m_fix_time')
row_update <- 1

######### MODEL 4 ##########
for (sf_val in sf_list) {
  for (sm_val in sm_list) {
    for (sz_val in sz_list) {
      # Specify initial genotype frequencies (4 female, 6 male)
      Z = 0.99
      z = 0.01
      M = 0.99
      m = 0.01
      
      current_z <- z
      current_m <- m
      
      z_fix <- FALSE
      m_fix <- FALSE
      
      #Initialize genotype frequencies
      #Females
      ZWM = Z * M
      zWM = z * M
      ZWm = Z * m
      zWm = z * m
      
      #Males
      ZZM = Z * Z * M
      ZzM = 2 * Z * z * M
      zzM = z * z * M
      ZZm = Z * Z * m
      Zzm = 2 * Z * z * m
      zzm = z * z * m
      
      for (gen in 1:generations) {
        if (current_z > 0.9999 & z_fix==FALSE) {
          z_fix_gen <- gen
          z_fix <- TRUE
        }
        if (current_m > 0.9999 & m_fix==FALSE) {
          m_fix_gen <- gen
          m_fix <- TRUE
        }
        if ((z_fix == TRUE) & (m_fix == TRUE)) {
          break
        }
        # Weight each male and female genotype by fitness and normalize
        Fsum= (ZWM + zWM) + (1+sf_val)*(ZWm + zWm)
        ZWM=ZWM/Fsum
        zWM=zWM/Fsum
        ZWm=(ZWm*(1+sf_val))/Fsum
        zWm=(zWm*(1+sf_val))/Fsum
        
        
        Msum = (ZZM + ZzM + zzM) + (1 - sm_val)*ZZm + (1 - sm_val + sz_val/2)*Zzm + (1-sm_val+sz_val)*zzm
        ZZM = ZZM/Msum
        ZzM = ZzM/Msum
        zzM = zzM/Msum
        ZZm = (ZZm * (1-sm_val))/Msum
        Zzm = (Zzm * (1-sm_val+(sz_val/2)))/Msum
        zzm = (zzm * (1-sm_val+sz_val))/Msum
        
        # Calculate the egg and sperm frequencies
        
        eZM=.5*ZWM
        eZm=.5*ZWm
        ezM=.5*zWM
        ezm=.5*zWm
        eWM=.5*(ZWM+zWM)
        eWm=.5*(ZWm+zWm)
        
        sZ=(ZZM + ZZm) +.5*(ZzM + Zzm)
        sz=.5*(ZzM+Zzm) + (zzM + zzm)
        
        # Calculate zygote frequencies the next generstion (zygotes are the same freq in both sexes)
        
        ZWM=eWM*sZ
        zWM=eWM*sz
        ZWm=eWm*sZ
        zWm=eWm*sz
        
        #print(paste('F zy Freq: ',ZWM + zWM + ZWm + zWm))
        
        ZZM = eZM * sZ
        ZzM = eZM * sz + ezM * sZ
        zzM = ezM * sz
        ZZm = eZm * sZ
        Zzm = eZm * sz + ezm * sZ
        zzm = ezm * sz
        
        current_z_female <- (zWM + zWm)/(ZWM + zWM + ZWm + zWm)
        current_z_male <- (ZzM + 2*zzM + Zzm + 2*zzm)/(2*(ZZM + ZzM + zzM + ZZm + Zzm + zzm))
        current_z <- (current_z_male + current_z_female)/2
        
        current_m_female <- (ZWm + zWm)/(ZWM + zWM + ZWm + zWm)
        current_m_male <- (ZZm + Zzm + zzm)/(ZZM + ZzM + zzM + ZZm + Zzm + zzm)
        current_m <- (current_m_female + current_m_male)/2
      }
      if (z_fix==FALSE) {
        z_fix_gen <- generations
      }
      if (m_fix==FALSE) {
        m_fix_gen <- generations
      }
      fix_time.df[row_update,] <- c(round(sf_val,digits=3),round(sm_val,digits=3),round(sz_val,digits=3),z_fix_gen,m_fix_gen)
      row_update <- row_update + 1
    }
  }
}

pdf('./model4_fix_time_z_heatmap.pdf')
for (i in round(sf_list,digits=3)) {
  fix_time_heatmap <- dcast(fix_time.df[fix_time.df$sf==i,c(2,3,4)],sm~sz, value.var = 'z_fix_time')
  row.names(fix_time_heatmap) <- fix_time_heatmap$sm
  fix_time_heatmap <- as.matrix(fix_time_heatmap[,-1])
  if (!all(fix_time_heatmap[1,1]==fix_time_heatmap)) {
    title <- paste("z Fix Time For sf = ", i)
    heatmap.2(fix_time_heatmap,Rowv=FALSE,Colv=FALSE,dendrogram='none',main=title,xlab='Z Chromosome',ylab='sm',key=T,trace='none',cexCol=1.25, cexRow=1.25,breaks=seq(0,1000,50),col=redblue(20))
  }
  else {
    print(paste("z Fix Time the Same for sf= ", i))
  }
}
dev.off()

pdf('./model4_fix_time_m_heatmap.pdf')
for (i in round(sf_list,digits=3)) {
  fix_time_heatmap <- dcast(fix_time.df[fix_time.df$sf==i,c(2,3,5)],sm~sz, value.var = 'm_fix_time')
  row.names(fix_time_heatmap) <- fix_time_heatmap$sm
  fix_time_heatmap <- as.matrix(fix_time_heatmap[,-1])
  if (!all(fix_time_heatmap[1,1]==fix_time_heatmap)) {
    title <- paste("m Fix Time For sf = ", i)
    heatmap.2(fix_time_heatmap,Rowv=FALSE,Colv=FALSE,dendrogram='none',main=title,xlab='Z Chromosome',ylab='sm',key=T,trace='none',cexCol=1.25, cexRow=1.25,breaks=seq(0,1000,50),col=redblue(20))
  }
  else {
    print(paste("m Fix Time the Same for sf= ", i))
  }
}
dev.off()