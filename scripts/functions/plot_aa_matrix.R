plot_aa_matrix<- function(aa_matrix,row_names,col_names,add_loess=T,plot_title="",add_mean_column=F){
  aa_matrix_sorted<- aa_matrix[intersect(row_names,rownames(aa_matrix)),intersect(col_names,colnames(aa_matrix)),drop=F]
  cellcolors<-matrix(NA,nrow=nrow(aa_matrix_sorted),ncol=ncol(aa_matrix_sorted))
  if(add_mean_column) aa_matrix_sorted<- cbind(aa_matrix_sorted,mean=rowMeans(aa_matrix_sorted))
  if(sum(aa_matrix_sorted<0)==0) cellcolors<- color.scale(aa_matrix_sorted,cs1=c(1,0),cs2=c(1,0),cs3=c(1,0))
  else{
    cellcolors[aa_matrix_sorted >= 0]<-
      color.scale(aa_matrix_sorted[aa_matrix_sorted >= 0],
                  cs1=c(1,1),cs2=c(1,0),cs3=c(1,0))
    cellcolors[aa_matrix_sorted < 0]<-
      color.scale(aa_matrix_sorted[aa_matrix_sorted < 0],
                  cs1=c(0,1),cs2=c(0,1),cs3=c(1,1))
  }
  color2D.matplot(aa_matrix_sorted,cellcolors=cellcolors,show.values=F,xlab =NA,ylab=NA,axes=FALSE,border=NA,main=plot_title)
  axis(2,at=0.5:(nrow(aa_matrix_sorted)-0.5),las=2,labels=rev(rownames(aa_matrix_sorted)),cex.axis=1)
  axis(1,labels = colnames(aa_matrix_sorted),at=(0.5:(ncol(aa_matrix_sorted)-0.5)),las=2,tick = FALSE)
  if(add_loess){
    if(sum(aa_matrix_sorted<0)==0){
      for(i in 1:nrow(aa_matrix_sorted)){
        y_coo_change<- nrow(aa_matrix_sorted)-i
        points(loess.smooth(1:ncol(aa_matrix_sorted),aa_matrix_sorted[i,])$x,y_coo_change+loess.smooth(1:ncol(aa_matrix_sorted),aa_matrix_sorted[i,])$y,col="red",type="l",ylim=c(y_coo_change,y_coo_change+1),lwd=3)
      }
    }
    else{
      for(i in 1:nrow(aa_matrix_sorted)){
        y_coo_change<- nrow(aa_matrix_sorted)-i+0.5
        points(loess.smooth(1:ncol(aa_matrix_sorted),aa_matrix_sorted[i,])$x,y_coo_change+loess.smooth(1:ncol(aa_matrix_sorted),aa_matrix_sorted[i,])$y,col="black",type="l",ylim=c(y_coo_change-1,y_coo_change+1),lwd=3)
      }
    }
  }  
  return(aa_matrix_sorted)
}
