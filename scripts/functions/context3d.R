## This script creates a "legoplot" similar to those produced by the Broad Institute
## The plot shows the relative abundance of each of the 6 possible mutations in the 
## 16 sequence contexts

## Load packages
library(rgl)

#### START OF FUNCTIONS

## Functions modified from the "demo(hist3d)" examples in the rgl package:
# library(rgl)
# demo(hist3d)
## Note; would it have killed the original author to comment their code?

## Draws a single "column" or "stack".
## X and Y coordinates determine the area of the column
## The Z coordinate determines the height of the column
## We include "lit=FALSE" arguments to remove the nasty shiny surfaces caused by lighting
stackplot.3d<-function(x,y,z,alpha=1,topcol="#078E53",sidecol="#aaaaaa"){
  
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to draw the sides and ends of the column separately  
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
  ## Determine the coordinates of each surface of the column and its edges
  x1=c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
  z1=c(rep(0,4),rep(c(0,0,z,z),4))
  y1=c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
  x2=c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2))
  z2=c(rep(c(0,z),4),rep(0,8),rep(z,8) )
  y2=c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2) )
  
  ## These lines create the sides of the column and its coloured top surface
  rgl.quads(x1,z1,y1,col=rep(sidecol,each=4),alpha=alpha,lit=FALSE)
  rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z,4),c(y[1],y[1],y[2],y[2]),
            col=rep(topcol,each=4),alpha=1,lit=FALSE) 
  ## This line adds black edges to the column
  rgl.lines(x2,z2,y2,col="#000000",lit=FALSE)
}
# Example:
# stackplot.3d(c(0,1),c(0,1),3,alpha=0.6)

## Calls stackplot.3d repeatedly to create a barplot
## z is the heights of the columns and must be an appropriately named vector
context3d<-function(z,alpha=1,scalexy=10,scalez=1,gap=0.2,plot_heatmap=F){
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to add each column sequentially
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
  ## Recreate Broad order
  types=c("C>G","T>A","C>A","T>G","C>T","T>C")
  # types=c("C>A","C>G","C>T","T>A","T>C","T>G")
  contexts=c("T_T","C_T","A_T","G_T","T_C","C_C","A_C","G_C",
             "T_A","C_A","A_A","G_A","T_G","C_G","A_G","G_G")
  contexts<- paste0(contexts,">",contexts)
  # typeorder=c()
  # for(type in types){
  #   typeorder=c(typeorder,paste(type,contexts,sep="_"))
  # }
  typeorder=c()
  for(type in types){
    contexts_tmp<- contexts
    substr(contexts_tmp,2,2)<- substr(type,1,1)
    substr(contexts_tmp,6,6)<- substr(type,3,3)
    typeorder=c(typeorder,contexts_tmp)
  }
  z=z[typeorder]
  
  ## Reorder data into 6 regions
  set1=c(1:4,17:20,5:8,21:24,9:12,25:28,13:16,29:32)
  set2=set1+32
  set3=set1+64
  neworder=c(set1,set2,set3)
  
  ## Define dimensions of the plot 
  dimensions=c(12,8)
  
  ## Scale column area and the gap between columns 
  y=seq(1,dimensions[1])*scalexy
  x=seq(1,dimensions[2])*scalexy
  gap=gap*scalexy
  
  ## Scale z coordinate
  z=z*scalez
  
  ## Set up colour palette
  # broadcolors=c("#805D3F","#72549A","#5EAFB2","#3F4F9D","#F2EC3C","#74B655")
  broadcolors=c("black","grey","blue","pink","red","green","red")
  colors=as.vector(sapply(broadcolors,rep,16))
  
  if(plot_heatmap){
    col_tmp<- colorRampPalette(colors = c("white","yellow","red"))(100)
    col_idx<- 50+log2(z/mean(z))*50
    col_idx[col_idx==0]<- 1
    col_idx[col_idx>100]<- 100
    colors=col_tmp[col_idx]
  }
  
  ## Plot each of the columns
  matrix_val<- matrix(NA,dimensions[2],dimensions[1],dimnames = list(rep(c("T","C","A","G"),dimensions[2]/4),rep(c("T","C","A","G"),dimensions[1]/4)))
  for(i in 1:dimensions[1]){
    for(j in 1:dimensions[2]){
      it=(i-1)*dimensions[2]+j # Variable to work out which column to plot; counts from 1:96
      stackplot.3d(c(gap+x[j],x[j]+scalexy),
                   c(-gap-y[i],-y[i]-scalexy),
                   z[neworder[it]],
                   alpha=alpha,
                   topcol=colors[neworder[it]],
                   sidecol=colors[neworder[it]])
      matrix_val[j,i]<- z[neworder[it]]
    }
  }
  ## Set the viewpoint and add axes and labels
  rgl.viewpoint(theta=50,phi=40,fov=0)
  axes3d("y-+",labels=TRUE)
  
  # return
  return(matrix_val)
}
# # Example:
# # context3d(counts)
# context3d(n_neg,scalexy = 1000,alpha=0.4)
# context3d(n_pos,scalexy = 1000,alpha=0.4)
# context3d(n_pos/n_neg,scalexy = 1000,scalez = 10000)
# 
# #### END OF FUNCTIONS
# 
# ## Read in example data and cast to an appropriate vector
# rawdata=read.table("snvspermegabase.txt",header=TRUE)
# counts=as.numeric(rawdata)
# names(counts)=colnames(rawdata)
# 
# ## Example plots
# context3d(counts)
# context3d(counts,alpha=0.4)

## Save your images to files if you wish
# rgl.snapshot(filename="example.png")

