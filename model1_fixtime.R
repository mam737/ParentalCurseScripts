# Model 1:  mito-auto Mother's curse - 11 Mar 2018
# Invasion and dynamics of Mother's curse mitochondrial type with autosomal restorer

# female genotypes AAM AaM aaM     AAm  Aam  aam
#   and fitnesses   1   1   1      1+sf 1+sf 1+sf

#   male genotypes AAMb AaMb aaMb     AAmb  Aamb        aamb
#   and fitnesses   1   1   1      1-sm 1-sm+sa/2  1-sm+sa

# four egg types:  AM aM   Am   am
# two sperm types: A   a

# Specify initial genotype frequencies (6 female, 6 male)

# Weight each male and female genotype by fitness and normalize

library(gplots)
library(reshape2)

## 3 paramters of interest, sf, sm, and sa
sf_list = seq(0,0.4,length.out = 20)
sm_list = seq(0,0.2,length.out = 20)
sa_list = seq(0,0.2,length.out = 20)
generations <- 1000

output.df <- data.frame(matrix(NA,nrow=(length(sf_list)*length(sm_list)*length(sa_list)), ncol=5))
colnames(output.df) <- c('sf','sm','sa' ,'a_fix_time','m_fix_time')
row_update <- 1

######### MODEL 1 ##########
for (sf_val in sf_list) {
  for (sm_val in sm_list) {
    for (sa_val in sa_list) {
      A = 0.99
      a = 0.01
      M = 0.99
      m = 0.01
      current_a <- a
      current_m <- m
      a_fix <- FALSE
      m_fix <- FALSE
      #Initialize genotype frequencies
      #Females
      AAM = A * A * M
      AAm = A * A * m
      AaM = 2 * A * a * M
      Aam = 2 * A * a * m
      aaM = a * a * M
      aam = a * a * m
      #Males
      AAMb = A * A * M
      AAmb = A * A * m
      AaMb = 2 * A * a * M
      Aamb = 2 * A * a * m
      aaMb = a * a * M
      aamb = a * a * m
      
      if (current_a != a | a_fix!=FALSE) {print("ERROR")}
      for (gen in 1:generations) {
        if (current_a > 0.9999 & a_fix==FALSE) {
          a_fix_gen <- gen
          a_fix <- TRUE
        }
        if (current_m > 0.9999 & m_fix==FALSE) {
          m_fix_gen <- gen
          m_fix <- TRUE
        }
        if ((a_fix == TRUE) & (m_fix == TRUE)) {
          break
        }
        # Weight each male and female genotype by fitness and normalize
        Fsum=(AAM+AaM+aaM)+(1+sf_val)*(AAm+Aam+aam)
        AAM=AAM/Fsum
        AaM=AaM/Fsum
        aaM=aaM/Fsum
        AAm=(AAm*(1+sf_val))/Fsum
        Aam=(Aam*(1+sf_val))/Fsum
        aam=(aam*(1+sf_val))/Fsum
        
        Msum=(AAMb+AaMb+aaMb)+(1-sm_val)*AAmb+(1-sm_val+sa_val/2)*Aamb+(1-sm_val+sa_val)*aamb
        AAMb=AAMb/Msum
        AaMb=AaMb/Msum
        aaMb=aaMb/Msum
        AAmb=(AAmb*(1-sm_val))/Msum
        Aamb=(Aamb*(1-sm_val+sa_val/2))/Msum
        aamb=(aamb*(1-sm_val+sa_val))/Msum
        
        # Calculate the egg and sperm frequencies
        
        eAM=AAM+.5*AaM
        eAm=AAm+.5*Aam
        eaM=aaM+.5*AaM
        eam=aam+.5*Aam
        
        sA=AAMb+AAmb+.5*AaMb+.5*Aamb
        sa=aaMb+aamb+.5*AaMb + .5*Aamb
        
        # Calculate zygote frequencies the next generstion (zygotes are the same freq in both sexes)
        
        AAM=eAM*sA
        AaM=eAM*sa+eaM*sA
        aaM=eaM*sa
        AAm=eAm*sA
        Aam=eAm*sa+eam*sA
        aam=eam*sa
        
        AAMb=eAM*sA
        AaMb=eAM*sa+eaM*sA
        aaMb=eaM*sa
        AAmb=eAm*sA
        Aamb=eAm*sa+eam*sA
        aamb=eam*sa
        
        current_a_female <- (2*aaM + AaM + 2*aam + Aam)/(2*(AAM + AaM + aaM + AAm + Aam + aam))
        current_a_male <- (2*aaMb + AaMb + 2*aamb + Aamb)/(2*(AAMb + AaMb + aaMb + AAmb + Aamb + aamb))
        current_a <- (current_a_female + current_a_male)/2
        
        current_m_female <- (AAm + Aam + aam)/(AAM + AaM + aaM + AAm + Aam + aam)
        current_m_male <- (AAmb + Aamb + aamb)/(AAMb + AaMb + aaMb + AAmb + Aamb + aamb)
        current_m <- (current_m_female + current_m_male)/2
        
      }
      if (a_fix==FALSE) {
        a_fix_gen <- generations
      }
      if (m_fix==FALSE) {
        m_fix_gen <- generations
      }
      output.df[row_update,] <- c(round(sf_val,digits=3),round(sm_val,digits=3),round(sa_val,digits=3),a_fix_gen,m_fix_gen)
      row_update <- row_update + 1
    }
  }
}

pdf('./model1_fix_time_a_heatmap.pdf')
for (i in round(sf_list,digits=3)) {
  fix_time_heatmap <- dcast(output.df[output.df$sf==i,c(2,3,4)],sm~sa, value.var = 'a_fix_time')
  row.names(fix_time_heatmap) <- fix_time_heatmap$sm
  fix_time_heatmap <- as.matrix(fix_time_heatmap[,-1])
  if (!all(fix_time_heatmap[1,1]==fix_time_heatmap)) {
    title <- paste("Fix Time For sf = ", i)
    heatmap.2(fix_time_heatmap,Rowv=FALSE,Colv=FALSE,dendrogram='none',main=title,xlab='Autosome',ylab='sm',key=T,trace='none',cexCol=1.25, cexRow=1.25,breaks=seq(0,1000,50),col=redblue(20))
  }
  else {
    print(paste("a Fix Time the Same for sf= ", i))
  }
}
dev.off()

pdf('./model1_fix_time_m_heatmap.pdf')
for (i in round(sf_list,digits=3)) {
  fix_time_heatmap <- dcast(output.df[output.df$sf==i,c(2,3,5)],sm~sa, value.var = 'm_fix_time')
  row.names(fix_time_heatmap) <- fix_time_heatmap$sm
  fix_time_heatmap <- as.matrix(fix_time_heatmap[,-1])
  if (!all(fix_time_heatmap[1,1]==fix_time_heatmap)) {
    title <- paste("m Fix Time For sf = ", i)
    heatmap.2(fix_time_heatmap,Rowv=FALSE,Colv=FALSE,dendrogram='none',main=title,xlab='Autosome',ylab='sm',key=T,trace='none',cexCol=1.25, cexRow=1.25,breaks=seq(0,1000,50),col=redblue(20))
  }
  else {
    print(paste("Fix Time the Same for sf= ", i))
  }
}
dev.off()