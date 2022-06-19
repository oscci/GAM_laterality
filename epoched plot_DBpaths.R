#---------------------------------------------------------------------------------------------------------------#
# GAM model (without contrasts) for LISA data (Word Generation) - simple stimulus (stim1 and stim2) - both runs
#
# Fixed single Gamma HRF

#For details of the general approach see Pedersen, E. J., Miller, D. L., Simpson, G. L., & Ross, N. (2019).
# Hierarchical generalized additive models in ecology: An introduction with mgcv. PeerJ, 7, e6876. https://doi.org/10.7717/peerj.6876

# All R scripts and datasets are available at the Open Science Framework repository: https://osf.io/gw4en/
  

#--------------------------------------------------------------------------------------------------#

# Created by Paul Thompson and Zoe Woodhead - 17th Oct 2019
# Edited by Paul Thompson - 26th Nov 2019
# Edited by Paul Thompson - 30th June 2020
# EDited by Paul Thompson - 16th July 2020
# EDited by Paul Thompson - 12th MAY 2021
# Edited by DVM Bishop - 19th June 2022 ; added comments; minor tweaks to data processing; 
#      modified GAM to include relative time/epoch
#--------------------------------------------------------------------------------------------------#

# install required packages for fitting the model.

#list_of_packages<-c("remotes","tidyverse","papaja","officer","fmri","knitr","utils","boot","ggpubr","psych","Rcpp","cladoRcpp","nlme","plm")
#new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
#if(length(new.packages))install.packages(new.packages,dependencies = TRUE)

#remotes::install_github("centerforopenscience/osfr")
#--------------------------------------------------------------------------------------------------#

## Packages
library(osfr)
library(utils)
require(dplyr)
require(tidyverse)
require(boot)
require(fmri)
require(ggpubr)
library(psych)
#library(nlme)
library(plm)
require(mgcv)
library(blandr)
library(gratia)
library(performance)
library(GGally)
library(here)
#--------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------#

order <- 3  #DB added; needed for creating gamma function

# This is the main function to run the analysis. The function does the following in order:
#
# PART 1:
#
#   Script takes a raw .exp datafile and preprocesses it ready for GAM analysis:
#    - Downsamples from 100 Hz to 25 Hz - 
#    - Identifies extreme outlying points that correspond to spikes or signal dropout and substitutes the mean for these values
#    - Identifies markers in the marker channel - these define epoch timings - 
#    - Creates a box car function showing when the task was ON or OFF - this is done separately for the period during word generation and word report. 
#   - Normalises the L and R fTCD signals to a mean of 100 by dividing by respective channel mean. This adjusts for any constant differences between L and R that may relate to angle of insonation. 
#   - Performs heart cycle integration (identified regular peaks in waveform and averages over the peak-to-peak interval). This removes a major, systematic source of variability that is of no interest. - It creates a channel corresponding to the epoch number - It performs baseline correction for each epoch separately for L and R by subtracting the mean value during the baseline period from the signal across the whole epoch. This ensures L and R are equated at the start of the epoch. - It saves the processed data into a .csv file
#   - Creates a column indicating the epoch number
#   - Creates a column showing relative time within the epoch
#   - Creates a box car function corresponding to the original POI (not used in current analysis)
#   - Saves the processed data into a .csv file
#
# PART 2:
#
#   - runs the gam
#   - saves the parameter estimates to data.frame.
#
# PART 3:
#
#   - plot the correlogram, Bland Altman plots, time series plots.
#

fTCD_gamma_LISA_GAM<-function(path,order)
{
  # get all files names to be loaded in and preprocessed
  filename1<-list.files(path,pattern = '.exp')
  
  ## Set parameters
  samplingrate <- 25 # Sampling rate after downsampling. Raw data is 100Hz, we take 1 in every 4 samples
  heartratemax <- 125
  
  # set up data.frame to hold the outputted parameter estimates from the GLMs.
  glm.data<-data.frame(matrix(NA,nrow=length(filename1),ncol=11))
  names(glm.data)<-c('ID',paste0('param',1:5),'HRF','AIC_glm','AIC_gam','BIC_glm','BIC_gam')
  
  #------------------------------------------------------------------------------------------------#
  ###########################################
  # PART 1: based on original fTCD analysis #
  #                                         #
  # Created by z.woodhead 30th July 2019    #
  # Edited  by z. woodhead 3rd Oct 2019     #
  ###########################################
  #------------------------------------------------------------------------------------------------#
  #Loop through to fit single GAM for each individual.
  for(j in 1:length(filename1))
 #   for(j in 1:10)
  {
    print(filename1[j])
    
    ## Read in raw data
    
    myfile <- filename1[j]
    mydata<-read.table(paste0(path,"/",myfile), skip = 6,  header =FALSE, sep ='\t')
    
    wantcols = c(2,3,4,7) #sec, L, R,marker #select columns of interest to put in shortdat
    shortdat = data.frame(mydata[,wantcols])
    rawdata = filter(shortdat, row_number() %% 4 == 0) # downsample to 5 Hz by taking every 20th point
    allpts = nrow(rawdata) # total N points in long file
    rawdata[,1] = (seq(from=1,to=allpts*4,by=4)-1)/100 #create 1st column which is time in seconds from start
    colnames(rawdata) = c("sec","L","R","marker")
    
    #----------------------------------------------------------
    ## Find markers; place where 'marker' column goes from low to high value
    # Marker channel shows some fluctuation but massive increase when marker is on so easy to detect
    
    mylen = nrow(rawdata); # Number of timepoints in filtered data (rawdata)
    markerplus = c(rawdata$marker[1] ,rawdata$marker); # create vectors with offset of one
    markerchan = c(rawdata$marker,0); 
    markersub = markerchan - markerplus; # start of marker indicated by large difference between consecutive data points
    meanmarker <- mean(rawdata$marker) # We will identify big changes in marker value that are > 5 sds
    markersize <- meanmarker+4*sd(rawdata$marker)
    origmarkerlist = which(markersub>markersize)
    norigmarkers = length(origmarkerlist)
    
    # Stimulus timings: Word Gen starts 5 seconds after marker and continues for 20 seconds (including REPORT phase)
    # Edit: stim1 models covert word generation, which starts 5 seconds after marker and continues for 15 seconds
    # stim2 models overt word reporting, which starts 20 seconds after marker and continues for 5 seconds
    stim1_delay_sec <- 5
    stim1_delay_samples <- stim1_delay_sec * samplingrate
    stim1_length_sec <- 15
    stim1_length_samples <- stim1_length_sec * samplingrate
    
    stim2_delay_sec <- 20
    stim2_delay_samples <- stim2_delay_sec * samplingrate
    stim2_length_sec <- 5
    stim2_length_samples <- stim2_length_sec * samplingrate
    
    rest_length_sec <- 30
    rest_length_samples <- rest_length_sec * samplingrate
    
    rawdata$stim1_on <- 0
    rawdata$stim2_on <- 0
    for (m in 1:norigmarkers){
      rawdata$stim1_on[(origmarkerlist[m]+stim1_delay_samples):(origmarkerlist[m]+stim1_delay_samples+stim1_length_samples)] <- 1
      rawdata$stim2_on[(origmarkerlist[m]+stim2_delay_samples):(origmarkerlist[m]+stim2_delay_samples+stim2_length_samples)] <- 1
    }
    #---------------------------------------------------------- (ZW NEW)
    # Identify raw datapoints below .0001 quartile (dropout_points) and above .9999 quartile (spike_points)
    
    dropout_points <- c(which(rawdata$L < quantile(rawdata$L, .0001)), 
                        which(rawdata$R < quantile(rawdata$R, .0001)))
    
    spike_points <- c(which(rawdata$L > quantile(rawdata$L, .9999)),
                      which(rawdata$R > quantile(rawdata$R, .9999)))
    
    #DB_June : in original script these points were identified but not smoothed/deleted prior to normalisation.
    #DB modified script to omit them from computation of mean (though it doesn't make much difference as they are so rare)
    
    #----------------------------------------------------------
    # Data normalisation
    
    meanL=mean(rawdata$L[-c(dropout_points,spike_points)])
    meanR=mean(rawdata$R[-c(dropout_points,spike_points)])
    rawdata$normal_L=rawdata$L/meanL * 100 
    rawdata$normal_R=rawdata$R/meanR * 100
    #For the dropout and spiking timepoints, substitute the mean (added by DB)
    rawdata$normal_L[c(dropout_points,spike_points)]<-meanL
    rawdata$normal_R[c(dropout_points,spike_points)]<-meanR
    
    
    
    #----------------------------------------------------------
    # Heartbeat integration
    peaklist=numeric(0)
    pdiff=numeric(0)
    badp=numeric(0)
    
    # Look through every sample from 6, to number of samples minus 6
    for(i in seq(6,mylen-6))
    {if(
      (rawdata$L[i] > rawdata$L[i-5])
      & (rawdata$L[i] > rawdata$L[i-4])
      & (rawdata$L[i] > rawdata$L[i-3])
      & (rawdata$L[i] > rawdata$L[i-2])
      & (rawdata$L[i] > rawdata$L[i-1])
      & (rawdata$L[i] > rawdata$L[i+1])
      & (rawdata$L[i] > rawdata$L[i+2])
      & (rawdata$L[i] > rawdata$L[i+3])
      & (rawdata$L[i]> rawdata$L[i+4])
      & (rawdata$L[i]> rawdata$L[i+5]))
    {peaklist=c(peaklist,i)
    }
    }
    
    # Check that the heartbeats are spaced by far enough!
    peakdiffmin = 60/heartratemax * samplingrate
    pdiff <- peaklist[2:length(peaklist)]-peaklist[1:(length(peaklist)-1)] # pdiff is a list of the number of samples between peaks
    badp<-which(pdiff<peakdiffmin) # badp is a list of the pdiff values that are less than peakdiffmin
    if (length(badp) != 0)
    {peaklist<-peaklist[-(badp+1)] # update peaklist, removing peaks identified by badp
    }
    #print(dim(rawdata))
    #print(peaklist)
    # Do heart beat integration
    peakn=length(peaklist)
    rawdata$heartbeatcorrected_L <- 0
    rawdata$heartbeatcorrected_R <- 0 
    for (p in 1:(peakn-1))
    {myrange=seq(peaklist[p],peaklist[p+1]) # the indices where the heartbeat will be replaced
    thisheart_L=mean(rawdata$normal_L[myrange]) # the new values that will be replaced
    thisheart_R=mean(rawdata$normal_R[myrange])
    rawdata$heartbeatcorrected_L[peaklist[p] : peaklist[p+1]]=thisheart_L
    rawdata$heartbeatcorrected_R[peaklist[p] : peaklist[p+1]]=thisheart_R
    if (p==1){
      rawdata$heartbeatcorrected_L[1:peaklist[p]] <- thisheart_L
      rawdata$heartbeatcorrected_R[1:peaklist[p]] <- thisheart_R
    }
    if (p==peakn-1){
      rawdata$heartbeatcorrected_L[peaklist[p] : mylen] <- thisheart_L
      rawdata$heartbeatcorrected_R[peaklist[p] : mylen] <- thisheart_R
    }
    }
    
    #-----------------------------------------------------------------------------------------------#
    # Save processed file in csv format
    #mynewfile <- paste0(getwd(),"/Lisa_data/Chpt4_fTCD_WordGen_rawdata/",strsplit(myfile, '*.exp'), '_processed.csv')
    #write.csv(rawdata, mynewfile, row.names=F)
    
    #-----------------------------------------------------------------------------------------------#
    #---------------------------------------------------------- (ZW NEW)
    # Identify extreme datapoints with values below 60 and above 140
    
    extreme_points <- c(which(rawdata$heartbeatcorrected_L < 60),
                        which(rawdata$heartbeatcorrected_L > 140),
                        which(rawdata$heartbeatcorrected_R < 60),
                        which(rawdata$heartbeatcorrected_R > 140))
    
    #DB new: create columns showing epoch and time relative to epoch start for each epoch (see below)
    rawdata$epoch<-NA #initialise new column
    rawdata$relativetime<-NA #initialise new column
    
    # Epoch timings
    epochstart_time   <- -12
    epochend_time     <- 40.52 #full epoch duration is 52.52 #updated by DB to include all of rest period
    epochstart_index  <- epochstart_time * samplingrate
    epochend_index    <- epochend_time * samplingrate
    basestart_time    <- -10 # baseline start
    baseend_time      <- 0 # baseline end
    basestart_index   <- basestart_time * samplingrate
    baseend_index    <- baseend_time * samplingrate
    
    # myepoched will be the full epoched trial
    myepoched <- array(0, dim=c(norigmarkers,epochend_index - epochstart_index + 1, 2)) # Set up an empty matrix
    
    for(mym in 1:norigmarkers) # for trials
    { 
      index1 = origmarkerlist[mym] + epochstart_index # index1 is index of the timepoint at the start of the epoch
      index2 = origmarkerlist[mym] + epochend_index # index2 is the index of the timepoint at the end of the epoch
      
      rawdata$relativetime[index1:index2]<-seq(from=epochstart_time, to=epochend_time, by=.04)
      
      
      # If recording started late, the start of the epoch for trial 1 will be beyond the recorded range. 
      # If this doesn't affect the baseline period (ie, results will be unaffected), then replace with mean
      
      if (index1 > 0){
        myepoched[mym,,1]=rawdata$heartbeatcorrected_L[index1:index2] #L side
        myepoched[mym,,2]=rawdata$heartbeatcorrected_R[index1:index2]
      }
    }
    
    # Baseline correction - (added to rawdata by DB. In fact, we don't use this in final analysis, but it can be useful to have it here, as it is the signal used when computing the LI from fTCD).
    basepoints=(basestart_index-epochstart_index):(baseend_index-epochstart_index) #all baseline points within epoch
    
    rawdata$baselinedL<-NA
    rawdata$baselinedR<-NA
    for (mym in 1:norigmarkers)
    {
      
      basemeanL=mean(myepoched[mym,basepoints,1]) #last dim is 3, which is HB corrected
      basemeanR=mean(myepoched[mym,basepoints,2])
      
      myepoched[mym,,1]=100+myepoched[mym,,1]-basemeanL #last dim 4 is HB and baseline
      myepoched[mym,,2]=100+myepoched[mym,,2]-basemeanR
      index1 = origmarkerlist[mym] + epochstart_index # index1 is index of the timepoint at the start of the epoch
      index2 = origmarkerlist[mym] + epochend_index # index2 is the index of the timepoint at the end of the epoch
      rawdata$baselinedL[index1:index2]<-myepoched[mym,,1]
      rawdata$baselinedR[index1:index2]<-myepoched[mym,,2]
      rawdata$epoch[index1:index2]<-mym
    }
    
    
    # Average over trials
    ntime <- dim(myepoched)[2]
    myepoched_average <- data.frame(
      "Lmean" <- rep(1, ntime),
      "Rmean" <- rep(1, ntime),
      "LRdiff" <- rep(1, ntime))
    
    myepoched_average$Lmean <- apply(myepoched[ , , 1], c(2), mean)
    myepoched_average$Rmean <- apply(myepoched[ , , 2], c(2), mean)
    myepoched_average$LRdiff <- myepoched_average$Lmean - myepoched_average$Rmean
    
    # # Plot myepoched_average
    
    
    myepoched_average$time<-seq(from=epochstart_time, to=epochend_time, by=.04)
    
    
    #########################################
    # PART 2                                #
    #                                       #
    # Created by P.Thompson 17th Oct 2019   #
    # Edited by P.Thompson 18th Oct 2019    #
    # Edited by P.Thompson 30th June 2020   #
    #########################################
    #-----------------------------------------------------------------------------------------------#
    rawdata2<-rawdata
    #-----------------------------------------------------------------------------------------------#
    rawdata$heartbeatcorrected_L[dropout_points]<-NA
    rawdata$heartbeatcorrected_R[dropout_points]<-NA
    
    rawdata$heartbeatcorrected_L[spike_points]<-NA
    rawdata$heartbeatcorrected_R[spike_points]<-NA 
    
    rawdata$heartbeatcorrected_L[extreme_points]<-NA
    rawdata$heartbeatcorrected_R[extreme_points]<-NA
    
    # Adapted 'fmri.stimulus' function from the R package 'fmri'. This is a condensed version that only gives option of the gamma HRF and convolves the HRF to the stimuli specified earlier in this script
    
    fmri.stimulus.PT2<- function(scans = dim(rawdata)[1],stim=stim, onsets = c(1,1+which(diff(rawdata$stim1_on)!=0)), durations = 375, TR = 1/25,scale=1)
    {
      onsets <- onsets * TR
      durations <- durations * TR
      onsets <- onsets * scale
      durations <- durations * scale
      scans <- scans * TR * scale
      TR <- TR/scale
      no <- length(onsets)
      durations <- rep(durations, no)
      
      stimulus<-stim
      
      .gammaHRF <- function(t, par = NULL) {
        th <- 0.242 * par[1]
        1/(th * factorial(3)) * (t/th)^3 * exp(-t/th)
      }
      
      
      par <- floor((durations[1]/28)*4)
      
      y <- .gammaHRF(0:(durations[1] * scale)/scale, par) 
      
      stimulus <-  convolve(stimulus, rev(y), type = "open")
      stimulus <- stimulus[unique((scale:scans)%/%(scale^2 * TR)) * scale^2 * TR]/(scale^2 * TR)
      stimulus <- stimulus - mean(stimulus)
      return(stimulus)
    }  
    
    #-----------------------------------------------------------------------------------------------#
    #-----------------------------------------------------------------------------------------------#
    
    # Create convolved stimulus function with HRF (applying the new fmri.stimulus.PT2 function above)
    
    gamma1 = fmri.stimulus.PT2(scans = dim(rawdata)[1], stim=rawdata$stim1_on, onsets = c(1,1+which(diff(rawdata$stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1)
    
    gamma2 = fmri.stimulus.PT2(scans = dim(rawdata)[1], stim=rawdata$stim2_on, onsets = c(1,1+which(diff(rawdata$stim2_on)!=0)), durations = stim2_length_samples, TR = 1,scale=1)
    
    #-----------------------------------------------------------------------------------------------# 
    # Binds all the stimuli into one matrix to be read into the fmri.design function. This converts the data into a design matrix and adds in the drift terms according to the order argument specified by the user.
    gamma = as.matrix(cbind(gamma1,gamma2))
    gamma = rbind(gamma,gamma)
    
    #-----------------------------------------------------------------------------------------------#
    
    # We create the design matrix and bind them together to give the same design matrix for each side (left and right), so that the main effect of side can be modelled appropriately.
    my_des<-fmri.design(gamma, order = order)
    
    # Add a dummy variable for side (signal). This is either 0 or 1 for left and right respectively.
    my_des<-cbind(my_des,rep(1:0,each=length(gamma1)))
    
    # Add a dummy variable for side (signal). This is either 0 or 1 for left and right respectively.
    my_des<-cbind(my_des,rep(1:0,each=length(gamma1)))
    
    # Add interaction variable for side (signal*stim1).
    my_des[,8]<-my_des[,8]*my_des[,1]
    
    # Use the design matrix to finish constructing the data for each GLM. 
    
    my_des[,7]<-dplyr::recode(my_des[,7], `0` = -1L, `1` = 1L)
    #-----------------------------------------------------------------------------------------------#
    
    mydata<-data.frame(y=c(rawdata$heartbeatcorrected_L,rawdata$heartbeatcorrected_R),stim1=my_des[,1],stim2=my_des[,2],t=my_des[,4],signal=as.factor(my_des[,7]),stim1_signal=my_des[,8])
    
    
    mydata<-data.frame(y=c(rawdata$heartbeatcorrected_L,rawdata$heartbeatcorrected_R),stim1=c(gamma1,gamma1),stim2=c(gamma2,gamma2),t=c(rawdata$sec,rawdata$sec),
                       relativetime=c(rawdata$relativetime,rawdata$relativetime),
                       epoch=as.factor(c(rawdata$epoch,rawdata$epoch)),signal=as.factor(rep(c(1,-1),each=length(gamma1))))
    
    levels(mydata$signal)<-c("Left","Right")
    #-----------------------------------------------------------------------------------------------#
    
    # filter out replicates in the dependent variable relating to the heartbeat correction (artificially induces autocorrelation if left in). The original pattern is sawtooth waveform after preprocessing, so sampling a single observation from the replicated observations reduced the computational load to estimate the model without affecting the fit as the autocorrelation is removed.
    #print(dim(mydata)[1]/2+1)
    #print(max(peaklist))
    
    #mydata_test<-mydata %>% group_by(y) %>% filter(t==min(t))
    mydata<-mydata[c(peaklist,peaklist+(dim(mydata)[1]/2+1)),] #peaklist duplicated for L and R
    
    #print(dim(mydata_test))
    #print(dim(mydata_test2))
    #print(length(unique(mydata$y)))
    #print(ts.plot(mydata$y))#;abline(v=mydata$y[peaklist]))
    #new1<-ggplot(mydata,aes(y=y,x=t))+geom_line()+geom_vline(xintercept=mydata$t[c(peaklist,peaklist+(dim(mydata)[1]+1))],colour="red")
    #new2<-ggplot(mydata,aes(y=y,x=t))+geom_line()+geom_vline(xintercept=mydata$t[c(peaklist)],colour="green")
    #print(new1)
    #print(new2)
    #print("1")
    # set optimisation parameters 
    glsControl(optimMethod = "L-BFGS-B",maxIter = 100)
    
    # fit gam model 

    myfit <- gam(y~s(t)+s(relativetime)+s(relativetime,by=epoch)+stim1+stim2+signal+stim1*signal,data=mydata)
    #print("2")
    myfit2 <- glm(y~stim1+stim2+t+I(t^2)+I(t^3)+signal+stim1*signal,data=mydata) #GLM comparison
    #print("3")
    #check fit of the model
    #print(appraise(myfit))
    
    print(paste0("AIC_gam=",AIC(myfit)))
    print(paste0("AIC_glm=",AIC(myfit2)))
    print(paste0("BIC_gam=",BIC(myfit)))
    print(paste0("BIC_glm=",BIC(myfit2)))
    
    #print("4")
    #print(summary(myfit))
    #-----------------------------------------------------------------------------------------------#
    
    # Extract the parameter estimates and record them for later use. Data stored in data.frame called 'glm.data'.
    glm.data[j,1] <- strsplit(basename(myfile),'[.]')[[1]][1]
    
    glm.data[j,7] <- "gamma"
    
    glm.data[j,2:6] <- anova(myfit)$'p.coeff'#parameter coefficients; last is interaction
    
    glm.data[j,8] <- AIC(myfit2) #print(paste0("AIC_gam=",AIC(myfit)))
    glm.data[j,9] <- AIC(myfit) #print(paste0("AIC_glm=",AIC(myfit2)))
    glm.data[j,10] <- BIC(myfit2)#print(paste0("BIC_gam=",BIC(myfit)))
    glm.data[j,11] <- BIC(myfit)#print(paste0("BIC_glm=",BIC(myfit2)))
    
    #-----------------------------------------------------------------------------------------------#
    #setup plot data (wrangling data to work with plot)
    
    mydata<-mydata%>%drop_na()
    
    myplotdat<-data.frame(y=mydata$y,x=mydata$t,fitted=predict(myfit),Signal=mydata$signal)#$gam, type = "response"
    
    #plot the time series for each individual per signal.
    g3<-ggplot(myplotdat,aes(y=y,x=x,colour=Signal))+geom_point(aes(colour=Signal),alpha=0.4)+geom_line(aes(y=fitted))+theme_bw()+theme(text=element_text(size=14))+ ylab('Normalised CBFV (cm/s)') + xlab('time(s)')
    
    # as we are fitting in a loop and printing to file, we need to use 'print' function with ggplot.
    print(g3)
    
    #================================================================================================#
    
    
    g4_epoch<-ggplot(myplotdat[1:200,],aes(y=y,x=x,colour=Signal))+geom_point(aes(colour=Signal),alpha=0.4)+geom_line(aes(y=fitted))+theme_bw()+theme(text=element_text(size=14))+ ylab('Normalised CBFV (cm/s)_epoch') + xlab('time(s)')
    
    print(g4_epoch)
    
    #output data
    glm_data<-glm.data
  }
  
  return(glm_data)    
}

################################# END OF FUNCTION ###################################################


#--------------------------------------------------------------------------------------------------#
# RUN FUNCTION FOR ALL PARTICIPANT DATA FILES
#--------------------------------------------------------------------------------------------------#

#Set the order


#pdf(file = '/Users/dorothybishop/Rprojects/GAM_laterality/Lisa_data/Fixed_HRF/gamma/HRF_signals_plots_LISA_WG_GAM_gamma_Jan2022.pdf', onefile = TRUE) #print plots to file.
my_results_LISA_WG_GAM_gamma_Jan2022<-fTCD_gamma_LISA_GAM(path='/Users/dorothybishop/Rprojects/GAM_laterality/Lisa_data/Chpt4_fTCD_WordGen_rawdata',order=order)
#dev.off()

#---------------------------------------------------------------------------------------------------#
write.csv(my_results_LISA_WG_GAM_gamma_Jan2022,'/Users/dorothybishop/Rprojects/GAM_laterality/PAULresults_LISA_WG_GAM_gamma.csv',row.names = FALSE)

my_results_LISA_WG_GAM_gamma_Jan2022<-my_results_LISA_WG_GAM_gamma_Jan2022[complete.cases(my_results_LISA_WG_GAM_gamma_Jan2022), ]
#---------------------------------------------------------------------------------------------------#


my_results_LISA_WG_GAM_gamma_ex_Jan2022 <- my_results_LISA_WG_GAM_gamma_Jan2022

#---------------------------------------------------------------------------------------------------#
# --------------------------------------------------------------------------------------------------#

#load LI based on old Doppler analysis method

old_res<-read.csv("/Users/dorothybishop/Rprojects/GAM_laterality/Lisa_data/WordGen_results.csv")

old_res<-old_res %>% rename(ID=Filename)

compare_results_Jan2022<-merge(my_results_LISA_WG_GAM_gamma_ex_Jan2022,old_res,by='ID',all.x = T)

fmri_data <- read.csv('/Users/dorothybishop/Rprojects/GAM_laterality/Lisa_data/Chapter5_fMRI_data.csv')

# Identify factors
factor_variables <- c('group_cat', 'group_lat', 'sex', 'hand_self_report', 'hand_QHP_cat', 'hand_EHI_cat', 'Old_fTCD_wg_cat', 'Old_fTCD_pptt_cat', 'fTCD_wg_cat', 'fTCD_pptt_cat')

for (i in 1:length(factor_variables))
{factor_ind <- str_which(colnames(fmri_data), paste0('^',factor_variables[i]))
fmri_data[,factor_ind] <- as.factor(fmri_data[,factor_ind])}

# Relabel group_cat and sex factors for clarity
# NB: group_cat 0=typical; 1=atypical
# sex 0=male; 1=female
levels(fmri_data$group_cat) <- c('T', 'A')
levels(fmri_data$sex) <- c('M', 'F')


fmri_data<-fmri_data[,c('ID','fMRI_diff_wg_frontal','fMRI_diff_wg_temporal','fMRI_diff_wg_MCA')]

big.summary<-merge(compare_results_Jan2022,fmri_data,by='ID')

wc<-which(colnames(big.summary)%in% c('param5','LI','fMRI_diff_wg_frontal','fMRI_diff_wg_temporal','fMRI_diff_wg_MCA'))
names(big.summary)[wc] <- c('GAM.LI','origLI','fmri.frontal','fmri.temp','fmri.MCA')


#output plots (correlagram and Bland Altman.)

#jpeg(file = '/Users/dorothybishop/Rprojects/GAM_laterality/Lisa_data/Fixed_HRF/gamma/HRF_signals_plots_LISA_WG_GAM_gamma_correlations.jpg')
psych::pairs.panels(big.summary[,c("origLI","GAM.LI","fmri.frontal","fmri.temp","fmri.MCA")])
#GGally::ggpairs(compare_results2[,c('LI (mean diff)','LI (GAM model-based)',"fMRI LI (Frontal)","fMRI LI (temporal)","fMRI LI (MCA)")])+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))
#dev.off()


#==================================================================#
ggally_cor_New<-
  function (data, mapping, ..., stars = TRUE, method = "pearson", 
            use = "complete.obs", display_grid = FALSE, digits = 3, title_args = list(...), 
            group_args = list(...), justify_labels = "right", align_percent = 0.5, 
            title = "Corr", alignPercent = warning("deprecated. Use `align_percent`"), 
            displayGrid = warning("deprecated. Use `display_grid`")) 
  {
    if (!missing(alignPercent)) {
      warning("`alignPercent` is deprecated. Please use `align_percent` if alignment still needs to be adjusted")
      align_percent <- alignPercent
    }
    if (!missing(displayGrid)) {
      warning("`displayGrid` is deprecated. Please use `display_grid`")
      display_grid <- displayGrid
    }
    na.rm <- if (missing(use)) {
      NA
    }
    else {
      (use %in% c("complete.obs", "pairwise.complete.obs", 
                  "na.or.complete"))
    }
    GGally::ggally_statistic(data = data, mapping = mapping, na.rm = na.rm, 
                             align_percent = align_percent, display_grid = display_grid, 
                             title_args = title_args, group_args = group_args, justify_labels = justify_labels, 
                             justify_text = "left", sep = if ("colour" %in% names(mapping)) 
                               ": "
                             else ":\n", title = title, text_fn = function(x, y) {
                               if (GGally:::is_date(x)) {
                                 x <- as.numeric(x)
                               }
                               if (GGally:::is_date(y)) {
                                 y <- as.numeric(y)
                               }
                               corObj <- stats::cor.test(x, y, method = method, 
                                                         use = use)
                               cor_est <- as.numeric(corObj$estimate)
                               cor_txt <- formatC(cor_est, digits = digits, format = "f")
                               if (isTRUE(stars)) {
                                 cor_txt <- str_c(cor_txt, GGally::signif_stars(corObj$p.value))
                                 cor_CI <- as.numeric(corObj$conf.int)
                                 cor_txt2 <- formatC(cor_CI, digits = digits, format = "f")
                                 cor_txt <- str_c(cor_txt, paste0("[",cor_txt2[1],',',cor_txt2[2],"]"),sep="\n")
                               }
                               cor_txt
                             })
  }

#==================================================================#


#psych::pairs.panels(compare_results2[,c('LI (mean diff)','LI (GAM model-based)',"fMRI LI (Frontal)","fMRI LI (temporal)","fMRI LI (MCA)")])
GGally::ggpairs(big.summary[,c('origLI','GAM.LI',"fmri.frontal","fmri.temp","fmri.MCA")], upper = list(continuous = ggally_cor_New))+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))+theme(text = element_text(size=14))

#ggsave(file = '/Users/dorothybishop/Rprojects/GAM_laterality/Lisa_data/Fixed_HRF/gamma/HRF_signals_plots_LISA_WG_GAM_gamma_correlations_ggally.jpg',width = 10,height = 10,dpi=600)

#jpeg(file = '/Users/dorothybishop/Rprojects/GAM_laterality/Lisa_data/Fixed_HRF/gamma/Bland_altman_plots_LISA_WG_GAM_gamma.jpg')
#print(blandr::blandr.draw(big.summary$'LI (new)', big.summary$'LI (old)') + theme_bw())
#dev.off()

# output demographics breakdown for summary in paper.
demographics<-read.csv("https://osf.io/x93w4/download",header=T)

demographics$ID<-substring(demographics$ID,1,7)

res.demographs<-merge(big.summary,demographics,by = "ID")

table(res.demographs$sex)

