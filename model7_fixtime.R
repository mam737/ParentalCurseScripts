# Model 7:  Y - X  Father's curse - 12 Mar 2018
# Invasion and dynamics of Father's curse X type with Y-linked restorer

# female genotypes  XX   Xx       xx     
#   and fitnesses   1   1-sf/2   1-sf      

#   male genotypes  XY     xY      Xy      xy
#   and fitnesses   1      1      1     1+sm  

# two egg types:     X and  x
# four sperm types:  X   x    Y    y


library(gplots)
library(reshape2)

## 2 paramters of interest, sf and sm
sf_list = seq(0,0.02,length.out = 20)
sm_list = seq(0,0.4,length.out = 20)
generations <- 2500

fix_time.df <- data.frame(matrix(NA,nrow=(length(sf_list)*length(sm_list)), ncol=4))
colnames(fix_time.df) <- c('sf','sm','y_fix_time','x_fix_time')
row_update <- 1

######### MODEL 7 ##########
for (sf_val in sf_list) {
  for (sm_val in sm_list) {
    # Specify initial genotype frequencies (3 female, 4 male)
    X = 0.99
    x = 0.01
    Y = 0.99
    y = 0.01
    current_x <- x
    current_y <- y  
    x_fix <- FALSE
    y_fix <- FALSE
    
    #Initialize genotype frequencies
    #Females
    XX = X * X
    Xx = 2 * X * x
    xx = x * x
    #Males
    XY = X * Y
    xY = x * Y
    Xy = X * y
    xy = x *y
    
    for (gen in 1:generations) {
      if (current_y > 0.9999 & y_fix==FALSE) {
        y_fix_gen <- gen
        y_fix <- TRUE
      }
      if (current_x > 0.9999 & x_fix==FALSE) {
        x_fix_gen <- gen
        x_fix <- TRUE
      }
      if ((y_fix == TRUE) & (x_fix == TRUE)) {
        break
      }
      # Weight each male and female genotype by fitness and normalize
      
      Fsum=XX+(1-sf_val/2)*Xx+(1-sf_val)*xx
      XX=XX/Fsum
      Xx=Xx*(1-sf_val/2)/Fsum
      xx=xx*(1-sf_val)/Fsum
      
      Msum=(XY+xY+Xy)+(1+sm_val)*xy
      XY=XY/Msum
      xY=xY/Msum
      Xy=Xy/Msum
      xy=(xy*(1+sm_val))/Msum
      
      # Calculate the egg and sperm frequencies
      
      eX=XX+.5*Xx
      ex=xx+.5*Xx
      
      sX=.5*(XY + Xy)
      sx=.5*(xY + xy)
      sY=.5*(XY + xY)
      sy=.5*(Xy + xy)
      
      # Calculate zygote frequencies the next generstion (zygotes are the same freq in both sexes)
      
      XX=eX*sX
      Xx=eX*sx+ex*sX
      xx=ex*sx

      
      XY=eX*sY
      xY=ex*sY 
      Xy=eX*sy
      xy=ex*sy
      
      current_x_female <- (2*xx + Xx)/(2*(XX+Xx+xx))
      current_x_male <- (xY + xy)/(XY+xY+Xy+xy)
      current_x <- (current_x_female + current_x_male)/2
      
      current_y <- (Xy + xy)/(XY + xY + Xy + xy)
    }
    if (y_fix==FALSE) {
      y_fix_gen <- generations
    }
    if (x_fix==FALSE) {
      x_fix_gen <- generations
    }
    fix_time.df[row_update,] <- c(round(sf_val,digits=3),round(sm_val,digits=3),y_fix_gen,x_fix_gen)
    row_update <- row_update + 1
  }
}

pdf('./model7_fix_time_y_heatmap.pdf')
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

pdf('./model7_fix_time_x_heatmap.pdf')
x_fix_time_heatmap <- dcast(fix_time.df[,c(1,2,4)],sf~sm, value.var ='x_fix_time')
row.names(x_fix_time_heatmap) <- x_fix_time_heatmap$sf
x_fix_time_heatmap <- as.matrix(x_fix_time_heatmap[,-1])
if (!all(x_fix_time_heatmap[1,1]==x_fix_time_heatmap)) {
  title = 'Fixation Time for x allele'
  heatmap.2(x_fix_time_heatmap,Rowv=FALSE, Colv=FALSE,dendrogram='none',main=title,xlab='Benefit to Males',ylab='Cost to Females',key=T,key.title='Fixation Time',key.xlab='Generations',density.info='none',trace='none',cexCol=1.25,cexRow=1.25,breaks=seq(0,2500,50),col=redblue(50))
} else {
  print(paste('Fixation Time for x allele is the Same Across all (sf,sm) values: ', x_fix_time_heatmap[1,1]))
}
dev.off()