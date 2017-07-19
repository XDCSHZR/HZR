# ---------------------- CNA simulation ----------------------------------------#

# ---------------------- parameters
# 1. Sample size (n)
# 2. Number of chromosomal sites or probes (m)
# 3. Number of recurrent CNAs (num_amp,num_del)
# 4. CN for each recurrent CNAs
# 5. Lengths of recurrent CNAs (len_amp,len_del)
# 6. Frequency of recurrent CNA gains and losses of the samples (freq_amp,freq_del)
# 7. Signal to noise ratio (SNR)

# Based on the above parameters, we simulate log2 intensity ratios


"simcna" <- function(n = 100, m = 2000)
{
 
for (loop in 1:10){
        # Initialization
        simu_data = matrix (0, n, m)
        simu_data[,]=2

        # Insert a number of recurrent CNA gains and losses at fixed positions across genome

        ap_s1=100
        ap_e1=149
        ap_s2=500
        ap_e2=529
        ap_s3=900
        ap_e3=919
        freq_amp = c(0.80,0.85,0.88)


           minsd = 0.6
           masd = 0.8


        # ----------------------------------------------- amplication ground truth
        num = freq_amp[1] * n
        for (i in 1:num){
             rs = round(runif(1, min = 1, max = n))         
             simu_data[rs,ap_s1:ap_e1] = 10
        }
  
        num = freq_amp[2] * n
        for (i in 1:num){
             rs = as.integer(runif(1, min = 1, max = n))
             simu_data[rs,ap_s2:ap_e2] = 4
        }

        num = freq_amp[3] * n
        for (i in 1:num){
             rs = as.integer(runif(1, min = 1, max = n))
             simu_data[rs,ap_s3:ap_e3] = 5
        }

        #-----------------------------------------------------------------------------
        
        #------ noise 1, random alteration
        for (i in 1:200){
             rand_s = as.integer(runif(1, min = 1, max = n-1))
             rand_pr = as.integer(runif(1, min = 1, max = m-500))
             CN = runif(1, min = 3, max = 4)
             len=round(runif(1,min=50,max=500))
             simu_data[rand_s,rand_pr:(rand_pr+len)] = CN
        }

        #------------------------------------------------- mixing of tumor cells and normal cells
        for (i in 1:n){
             rand_p = runif(1, min = 0.3, max = 0.7)
             simu_data[i,] = log2((simu_data[i,] * rand_p + 2*(1-rand_p))/2)
        }

        # ------------------------------------------------ add noise to the data
        for (i in 1 : n){
            rand_sd = runif(1, min = minsd, max = masd)      # generate random numbers
            Gauss_noise = rnorm(m, mean = 0, sd = rand_sd)
            simu_data[i,] = simu_data[i,] + Gauss_noise
        }


        # ------------------------------------------------ write out
        simu_data = t(simu_data)

        #-------------------------------------------------- used for cmds.
        header = matrix(0,1,n)
        header[1,] = 1:n
        simu_data_two=(2^simu_data)*2;
        simu_data_two = rbind(header, simu_data_two)
        header = matrix(0,m+1,2)
        header[1,1]='chromosome'
        header[1,2]='position'
        header[2:(m+1),1]=1
        header[2:(m+1),2]=1:m
        simu_data_two = cbind(header, simu_data_two)
        out.file_cmds = paste("D:\\RProject\\data\\simu_data",loop,sep="_")
        write.table(simu_data_two, file = out.file_cmds, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)   
   
        cat("the No.",loop, "has been simulated!","\n")  
    }
}

a<-function()
{
b=rnorm(100,4,0.4);
b=as.matrix(b);
out.file_cmds = paste("D:\\Teaching\\2015_course\\data_analysis\\sim_data\\1111",sep="_")
        write.table(b, file = out.file_cmds, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)   
}









