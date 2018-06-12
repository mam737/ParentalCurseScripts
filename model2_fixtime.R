# Model 1:  mito-auto Mother's curse - 11 Mar 2018
# Invasion and dynamics of Mother's curse mitochondrial type with autosomal restorer

# female genotypes AAM AaM aaM     AAm  Aam  aam
#   and fitnesses   1   1   1      1+sf 1+sf 1+sf

#   male genotypes AAMb AaMb aaMb     AAmb  Aamb        aamb
#   and fitnesses   1   1   1      1-sm 1-sm+sa/2  1-sm+sa

# four egg types:  AM aM   Am   am
# two sperm types: A   a

library(gplots)
library(reshape2)

## 3 paramters of interest, sf, sm, and sa
sf_list = seq(0,0.4,length.out = 20)
sm_list = seq(0,0.2,length.out = 20)
sx_list = seq(0,0.2,length.out = 20)
generations <- 1000

fix_time.df <- data.frame(matrix(NA,nrow=(length(sf_list)*length(sm_list)*length(sx_list)), ncol=5))
colnames(fix_time.df) <- c('sf','sm','sx' ,'x_fix_time','m_fix_time')
row_update <- 1

######### MODEL 1 ##########
for (sf_val in sf_list) {
  for (sm_val in sm_list) {
    for (sx_val in sx_list) {
      # Specify initial genotype frequencies (6 female, 6 male)
      X = 0.99
      x = 0.01
      M = 0.99
      m = 0.01
      current_x <- x
      current_m <- m
      x_fix <- FALSE
      m_fix <- FALSE
      #Initialize genotype frequencies
      #Females
      XXM = X * X * M
      XXm = X * X * m
      XxM = 2 * X * x * M
      Xxm = 2 * X * x * m
      xxM = x * x * M
      xxm = x * x * m
      #Males
      XM = X * M
      xM = x * M
      Xm = X * m
      xm = x * m
      
      for (gen in 1:generations) {
        if (current_x > 0.9999 & x_fix==FALSE) {
          x_fix_gen <- gen
          x_fix <- TRUE
        }
        if (current_x > 0.9999 & x_fix==FALSE) {
          x_fix_gen <- gen
          x_fix <- TRUE
        }
        if ((x_fix == TRUE) & (x_fix == TRUE)) {
          break
        }
        # Weight each male and female genotype by fitness and normalize
        Fsum=(XXM+XxM+xxM)+(1+sf_val)*(XXm+Xxm+xxm)
        XXM=XXM/Fsum
        XxM=XxM/Fsum
        xxM=xxM/Fsum
        XXm=(XXm*(1+sf_val))/Fsum
        Xxm=(Xxm*(1+sf_val))/Fsum
        xxm=(xxm*(1+sf_val))/Fsum
        
        Msum=(XM+xM)+(1-sm_val)*(Xm)+(1-sm_val+sx_val)*(xm)
        XM=XM/Msum
        xM=xM/Msum
        Xm=(Xm*(1-sm_val))/Msum
        xm=(xm*(1-sm_val+sx_val))/Msum
        
        # Calculate the egg and sperm frequencies
        
        eXM=XXM+.5*XxM
        eXm=XXm+.5*Xxm
        exM=xxM+.5*XxM
        exm=xxm+.5*Xxm
        
        sX=XM+Xm
        sx=xM+xm
        
        # Calculate zygote frequencies the next generstion (zygotes are the same freq in both sexes)
        
        XXM=eXM*sX
        XxM=eXM*sx+exM*sX
        xxM=exM*sx
        XXm=eXm*sX
        Xxm=eXm*sx+exm*sX
        xxm=exm*sx
        
        #print(paste('F zy Freq: ', XXM + XxM + xxM + XXm + Xxm + xxm))
        XM=eXM
        Xm=eXm
        xM=exM
        xm=exm
        
        current_x_female <- (2*xxM + XxM + 2*xxm + Xxm)/(2*(XXM + XxM + xxM + XXm + Xxm + xxm))
        current_x_male <- (xM + xm)/(XM + Xm + xM + xm)
        current_x <- (current_x_female + current_x_male)/2
        
        current_m_female <- (XXm + Xxm + xxm)/(XXM + XxM + xxM + XXm + Xxm + xxm)
        current_m_male <- (Xm + xm)/(XM + Xm + xM + xm)
        current_m <- (current_m_female + current_m_male)/2
      }
      if (x_fix==FALSE) {
        x_fix_gen <- generations
      }
      if (m_fix==FALSE) {
        m_fix_gen <- generations
      }
      fix_time.df[row_update,] <- c(round(sf_val,digits=3),round(sm_val,digits=3),round(sx_val,digits=3),x_fix_gen,m_fix_gen)
      row_update <- row_update + 1
    }
  }
}

pdf('./model2_fix_time_x_heatmap.pdf')
for (i in round(sf_list,digits=3)) {
  fix_time_heatmap <- dcast(fix_time.df[fix_time.df$sf==i,c(2,3,4)],sm~sx, value.var = 'x_fix_time')
  row.names(fix_time_heatmap) <- fix_time_heatmap$sm
  fix_time_heatmap <- as.matrix(fix_time_heatmap[,-1])
  if (!all(fix_time_heatmap[1,1]==fix_time_heatmap)) {
    title <- paste("x Fix Time For sf = ", i)
    heatmap.2(fix_time_heatmap,Rowv=FALSE,Colv=FALSE,dendrogram='none',main=title,xlab='X Chromosome',ylab='sm',key=T,trace='none',cexCol=1.25, cexRow=1.25,breaks=seq(0,1000,50),col=redblue(20))
  }
  else {
    print(paste("x Fix Time the Same for sf= ", i))
  }
}
dev.off()

pdf('./model2_fix_time_m_heatmap.pdf')
for (i in round(sf_list,digits=3)) {
  fix_time_heatmap <- dcast(fix_time.df[fix_time.df$sf==i,c(2,3,5)],sm~sx, value.var = 'm_fix_time')
  row.names(fix_time_heatmap) <- fix_time_heatmap$sm
  fix_time_heatmap <- as.matrix(fix_time_heatmap[,-1])
  if (!all(fix_time_heatmap[1,1]==fix_time_heatmap)) {
    title <- paste("m Fix Time For sf = ", i)
    heatmap.2(fix_time_heatmap,Rowv=FALSE,Colv=FALSE,dendrogram='none',main=title,xlab='X Chromosome',ylab='sm',key=T,trace='none',cexCol=1.25, cexRow=1.25,breaks=seq(0,1000,50),col=redblue(20))
  }
  else {
    print(paste("m Fix Time the Same for sf= ", i))
  }
}
dev.off()