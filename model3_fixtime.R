# Model 3: mito-Y Mother's curse - 11 March 2018
# Invasion and dynamics of Mother's curse mitochondrial type with Y-linked restorer

# female genotypes   M       m
#   and fitnesses    1       1+sf

#   male genotypes  YM      yM      Ym      ym
#   and fitnesses   1       1       1-sm    1-sm+sy

# two egg types:   M      m
# two sperm types: Y  and y

# Specify initial genotype frequencies (2 female, 4 male)

library(gplots)
library(reshape2)

## 3 paramters of interest, sf, sm, and sa
sf_list = seq(0,0.4,length.out = 20)
sm_list = seq(0,0.2,length.out = 20)
sy_list = seq(0,0.2,length.out = 20)
generations <- 1000

fix_time.df <- data.frame(matrix(NA,nrow=(length(sf_list)*length(sm_list)*length(sy_list)), ncol=5))
colnames(fix_time.df) <- c('sf','sm','sy' ,'y_fix_time','m_fix_time')
row_update <- 1

######### MODEL 3 ##########
for (sf_val in sf_list) {
  for (sm_val in sm_list) {
    for (sy_val in sy_list) {
      # Specify initial genotype frequencies (2 female, 4 male)
      Y = 0.99
      y = 0.01
      M_allele = 0.99
      m_allele = 0.01
      
      current_y <- y
      current_m <- m_allele
      y_fix <- FALSE
      m_fix <- FALSE
      
      #Initialize genotype frequencies
      #Females
      M = M_allele
      m = m_allele
      
      YM = Y * M_allele
      Ym = Y * m_allele
      yM = y * M_allele
      ym = y * m_allele
      
      for (gen in 1:generations) {
        if (current_y > 0.9999 & y_fix==FALSE) {
          y_fix_gen <- gen
          y_fix <- TRUE
        }
        if (current_m > 0.9999 & m_fix==FALSE) {
          m_fix_gen <- gen
          m_fix <- TRUE
        }
        if ((y_fix == TRUE) & (m_fix == TRUE)) {
          break
        }
        # Weight each male and female genotype by fitness and normalize
        Fsum=(M)+(1+sf_val)*(m)
        M=M/Fsum
        m=(m*(1+sf_val))/Fsum
        
        Msum=(YM+yM)+(1-sm_val)*(Ym)+(1-sm_val+sy_val)*(ym)
        YM=YM/Msum
        yM=yM/Msum
        Ym=(Ym*(1-sm_val))/Msum
        ym=(ym*(1-sm_val+sy_val))/Msum
        
        # Calculate the egg and sperm frequencies
        
        eM=M
        em=m
        
        sY=YM+Ym
        sy=yM+ym
        
        # Calculate zygote frequencies the next generstion (zygotes are the same freq in both sexes)
        
        M=eM
        m=em
        
        YM=eM*sY
        Ym=em*sY
        yM=eM*sy
        ym=em*sy
        
        current_y <- (yM + ym)/(YM+Ym+yM+ym)
        
        current_m_female <- (m)/(M+m)
        current_m_male <- (Ym + ym)/(YM + Ym + yM + ym)
        current_m <- (current_m_female + current_m_male)/2
      }
      if (y_fix==FALSE) {
        y_fix_gen <- generations
      }
      if (m_fix==FALSE) {
        m_fix_gen <- generations
      }
      fix_time.df[row_update,] <- c(round(sf_val,digits=3),round(sm_val,digits=3),round(sy_val,digits=3),y_fix_gen,m_fix_gen)
      row_update <- row_update + 1
    }
  }
}

pdf('./model3_fix_time_y_heatmap.pdf')
for (i in round(sf_list,digits=3)) {
  fix_time_heatmap <- dcast(fix_time.df[fix_time.df$sf==i,c(2,3,4)],sm~sy, value.var = 'y_fix_time')
  row.names(fix_time_heatmap) <- fix_time_heatmap$sm
  fix_time_heatmap <- as.matrix(fix_time_heatmap[,-1])
  if (!all(fix_time_heatmap[1,1]==fix_time_heatmap)) {
    title <- paste("y Fix Time For sf = ", i)
    heatmap.2(fix_time_heatmap,Rowv=FALSE,Colv=FALSE,dendrogram='none',main=title,xlab='Y Chromosome',ylab='sm',key=T,trace='none',cexCol=1.25, cexRow=1.25,breaks=seq(0,1000,50),col=redblue(20))
  }
  else {
    print(paste("y Fix Time the Same for sf= ", i))
  }
}
dev.off()

pdf('./model3_fix_time_m_heatmap.pdf')
for (i in round(sf_list,digits=3)) {
  fix_time_heatmap <- dcast(fix_time.df[fix_time.df$sf==i,c(2,3,5)],sm~sy, value.var = 'm_fix_time')
  row.names(fix_time_heatmap) <- fix_time_heatmap$sm
  fix_time_heatmap <- as.matrix(fix_time_heatmap[,-1])
  if (!all(fix_time_heatmap[1,1]==fix_time_heatmap)) {
    title <- paste("m Fix Time For sf = ", i)
    heatmap.2(fix_time_heatmap,Rowv=FALSE,Colv=FALSE,dendrogram='none',main=title,xlab='Y Chromosome',ylab='sm',key=T,trace='none',cexCol=1.25, cexRow=1.25,breaks=seq(0,1000,50),col=redblue(20))
  }
  else {
    print(paste("m Fix Time the Same for sf= ", i))
  }
}
dev.off()