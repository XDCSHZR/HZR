#source("http://bioconductor.org/biocLite.R")
#biocLite("DNAcopy")

library(DNAcopy)     # gistic,SCEICG,stac

"CBS_data"<-function(){
   data.dir="D:\\RProject\\data\\"
   ind=c(1:10)
   for(i in ind){
         cat("The ",i, "th data", "\n")
         input.data.dir = paste(data.dir,"simu_data_",i,sep="")
         data=read.table(input.data.dir)
         data=data[-1,-(1:2)]
         data=log2(data/2)
         head=matrix(0,nrow(data),3)
         head[,1]=1
         head[,2]=1:nrow(data)
         head[,3]=1:nrow(data)

         chrom <- rep(1,nrow(data))
         maploc <- 1:nrow(data)          
         seg.file_g=matrix(0,1,6)
         seg.file_g_one=matrix(0,1,6)
         seg.file=matrix(0,nrow(data),1)

         #---------------------------------------------------------------------------stac
         stac_amp = matrix(0,1,nrow(data))
         stac_amp[1,]=1:nrow(data)
         stac_amp_one=matrix(0,1,nrow(data))

         stac_del = matrix(0,1,nrow(data))
         stac_del[1,]=1:nrow(data)
         stac_del_one=matrix(0,1,nrow(data))

         for (j in 1:ncol(data)){
             cat("sampl No,",j,"\n")
             seg<- segment(CNA(data[,j],chrom,maploc))

             n=0
             for (k in 1:length(seg$output$loc.start)){
                  seg.file_g_one[1,1]=j
                  seg.file_g_one[1,2]=1
                  seg.file_g_one[1,3]=seg$output$loc.start[k]*10000
                  seg.file_g_one[1,4]=seg$output$loc.end[k]*10000
                  seg.file_g_one[1,5]=seg$output$num.mark[k]
                  seg.file_g_one[1,6]=seg$output$seg.mean[k]
                  seg.file_g=rbind(seg.file_g,seg.file_g_one)
                  seg.file_g_one=matrix(0,1,6)

                  len=seg$output$num.mark[k]
                  for(l in 1:len){
                      n = n+1
                      seg.file[n,1]=seg$output$seg.mean[k]
            
                      #------------------------------------------stac
                      if (as.numeric(seg$output$seg.mean[k])>0.1)
                          stac_amp_one[1,n]=1
                      else
                         stac_amp_one[1,n]=0

                      if (as.numeric(seg$output$seg.mean[k]) < (-0.1))
                          stac_del_one[1,n]=1
                      else
                          stac_del_one[1,n]=0
                  }
             }
             head = cbind(head,seg.file)

             stac_amp=rbind(stac_amp,stac_amp_one)
             stac_del=rbind(stac_del,stac_del_one)


         }
         seg.file_g = seg.file_g[-1,]
         out.file=paste(data.dir,"seg",i,sep="")
         write.table(seg.file_g,file=out.file,row.names=F,col.names=F,quote=F,sep="\t")

         out.file=paste(data.dir,"data_denoised_",i,sep="")
         write.table(head,file=out.file,row.names=F,col.names=F,quote=F,sep="\t")

         #-----------------------------------------------------------------------------------stac
         
         firstLine=stac_amp[1,]
         firstLine[1]=paste("\t",firstLine[1],sep="")
         firstLine=as.matrix(firstLine)
         firstLine=t(firstLine)
         stac_amp=stac_amp[-1,]
         stac_del=stac_del[-1,]

         h=matrix(0,nrow(stac_amp),1)
         h[,1]=1:nrow(stac_amp)
         stac_amp=cbind(h,stac_amp)
         stac_del=cbind(h,stac_del)


         out.file=paste(data.dir,"\\","data_stac","\\","data_stac_amp",i,sep="")
         #write.table(firstLine,file=out.file,row.names=F,col.names=F,quote=F,sep="\t")
         #write.table(stac_amp,file=out.file,append=TRUE,row.names=F,col.names=F,quote=F,sep="\t")

         #out.file=paste(data.dir,"\\","data_stac","\\","data_stac_del",i,sep="")
         #write.table(firstLine,file=out.file,row.names=F,col.names=F,quote=F,sep="\t")
         #write.table(stac_del,file=out.file,append=TRUE,row.names=F,col.names=F,quote=F,sep="\t")
   }
}

























