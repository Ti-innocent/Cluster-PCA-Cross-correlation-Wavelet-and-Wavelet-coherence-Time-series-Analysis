###################-------------------###################
#############                                 ###########
#               Preprosessing the Data                  #

#############                                 ###########
###################-------------------###################

#rm(list=ls())                         # this command clears the history

data<-read.csv('event2.csv')   # importing the data
View(data)
rownames(data) <- data[,1]            # Setting rownames

new_data=data[,-1]     # removing the first column because it only serves as our row names

class(new_data)        # checking your data class

dim(new_data)          # getting the dimension of your data
names(new_data)        # getting the column names

str(new_data)          # checking data structure

anyNA(new_data)        # checking for missing values

summary(new_data)      # summary of the data

View(new_data)         # view your data

###################-------------------###################
#############                                 ###########
#                    Cluster Analysis                   #
#############                                 ###########
###################-------------------###################


## standardizing and preparing the dataset
new_data_scaled <- scale(new_data)                # standardizing the variables
new_data_scaled.dist <- dist(new_data_scaled)     # distance computation with year as observation
new_data_scaled.distT <- dist(t(new_data_scaled)) # distance computation with year as variables

# performing hierarchical clustering of the observations (variables) using Ward method,
#single and average linkage 
ward = hclust(new_data_scaled.distT,method='ward.D2')

## Plotting the dendogram
#png('cluster_dendogram_variable.png')    # to put the figure in a file
par(mfrow = c(1 ,1))
plot(ward, main = 'Wards Method',xlab =" ", sub ="")
rect.hclust(ward, k =3, border = 'red') ## selecting three clusters

###################---------------------------###################
#################                                 ###############
#                 Principal Componenet Analysis                 #
#################                                 ###############
###################---------------------------###################

library(psych)

## Rotated PCA after extracting four components
pc2 <- principal(new_data_scaled, nfactors = 3, rotate = "varimax") 
scores<- as.data.frame(pc2$scores)  #scores after rotation
plot(scores$RC1,type ="o", xlab = 'year', ylab='RC1', main='Rotated Component 1',col="red")
#dev.off()


#png('pca_scores_rotated2.png')
plot(scores$RC2~rownames(scores),type ="o", xlab = 'year', ylab='RC2', main='Rotated Component 2',col="red")
#dev.off()

plot(scores$RC3~rownames(scores),type ="o", xlab = 'year', ylab='RC3', main='Rotated Component 3',col="red")

################## Additional codes


###################---------------------------###################
#################                                 ###############
#                       Time Series Analysis                    #
#              Autocorrelation and Cross correlation            #
#################                                 ###############
###################---------------------------###################
library(forecast)
library(TSA)
library(biwavelet)


needed_data <-scores[, c(1,2,3)] # choosing RC1,RC2 and RC3
View(needed_data)

### converting data to time series
#png('timeseries_plot.png')
prep <- ts(needed_data$RC1, start = c(2000), end = c(2020), frequency = 1)
prep1 <- ts(needed_data$RC2, start = c(2000), end = c(2020), frequency = 1)
prep2 <- ts(needed_data$RC3, start = c(2000), end = c(2020), frequency = 1)
#dev.off()

#png('Auto_correlation.png')
par(mfcol = c(1,1))
plot(prep, ylab = 'Yearly hazards', main = 'Time Series plot of RC1', col='black')
abline(reg = lm(prep ~ time(prep)))
plot(prep1, ylab = 'Yearly hazards', main = 'Time Series plot of RC2', col='black')
abline(reg = lm(prep1 ~ time(prep1)))
plot(prep2, ylab = 'Yearly hazards', main = 'Time Series plot of RC3', col='black')
abline(reg = lm(prep2 ~ time(prep2)))
#dev.off()
################################### detrend data and plot
tr1=lm(prep~c(1:length(prep)))
RC1=residuals(tr1)
tr2=lm(prep1~c(1:length(prep1)))
RC2=residuals(tr2)
tr3=lm(prep2~c(1:length(prep2)))
RC3=residuals(tr3)
length(RC3)
##### Plot detrended data set ##########
par(mfcol = c(1,1))
plot.ts(RC1, ylab = 'Yearly hazards', main = 'Time Series plot for RC1', col='black')
abline(reg = lm(prep ~ time(prep)))
plot.ts(RC2, ylab = 'Yearly hazards', main = 'Time Series plot for RC2', col='black')
plot.ts(RC3, ylab = 'Yearly hazards', main = 'Time Series plot for RC3', col='black')
#dev.off()

### Autocorrelation plot
#png('Auto_correlation.png')
par(mfcol = c(1,1))
forecast::Acf(RC1, lag.max = 20)
forecast::Acf(RC2, lag.max = 20)
forecast::Acf(RC3, lag.max = 20)
#dev.off()


### cross-correlation plot
#png('cross_correlation.png')
par(mfcol = c(1,1))
forecast::Ccf(RC1,RC2,lag.max =20)
forecast::Ccf(RC1,RC3,lag.max = 20)
forecast::Ccf(RC2,RC3,lag.max = 20)
#dev.off()

###################---------------------------###################
#################                                 ###############
#                       Time Series Analysis                    #
#                         Fourier Analysis                      #
#################                                 ###############
###################---------------------------###################
#png('periodogram.png')
par(mfcol = c(1,1))
p1 <- periodogram(RC1)
p2 <- periodogram(RC2)
p3 <- periodogram(RC3)
## converting the frequency to periods
#png('frequency_gram.png')
period_RC1 <- 1/p1$freq
period_RC2 <- 1/p2$freq
period_RC3 <- 1/p3$freq
#par(mfcol = c(1,3))

plot(period_RC1, p1$spec, type = 'b', xlim = c(1,21), xlab = 'period', ylab = 'periodogram',
     main = 'RC1')
plot(period_RC2, p2$spec, type = 'b', xlim = c(1,21), xlab = 'period', ylab = 'periodogram',
     main = 'RC2')
plot(period_RC3, p3$spec, type = 'b', xlim = c(1,21), xlab = 'period', ylab = 'periodogram',
     main = 'RC3')
#dev.off()
#################################################################
#                                                               #
#                       Time Series Analysis                    #
#             Wavelet Analysis and Wavelet Coherence            #
#################                                 ###############
#################################################################
## Wavelet
timelong=seq(2000, 2021-1, 1)  # indexing
wavelet1=wt(cbind(timelong,RC1),dj=0.1,mother="morlet",max.scale=8) #wavelet
wavelet2=wt(cbind(timelong,RC2),dj=0.1,mother="morlet",max.scale=8) #wavelet
wavelet3=wt(cbind(timelong,RC3),dj=0.1,mother="morlet",max.scale=8) #wavelet
## plot of wavelet
#png('wavelet.png')
#par(mfcol = c(1,1))
par(oma=c(0, 0, 0, 1), mar=c(5, 4, 4, 5) + 0.1)
xx=plot(wavelet1,  plot.cb = TRUE,lwd.coi = 1,col.coi="white",alpha.coi=0.5,col.sig="black",
        lwd.sig = 2, ncol = 768, tol = 0.95, plot.phase=F,ylab="Period(Year)", xlab = "Time(Year)",
        main='RC1')

#par(oma=c(0, 0, 0, 1), mar=c(5, 4, 4, 5) + 0.1)
yy=plot(wavelet2,  plot.cb = TRUE,lwd.coi = 1,col.coi="white",alpha.coi=0.5,col.sig="black",
        lwd.sig = 2, ncol = 768, tol = 0.95, plot.phase=F,ylab="Period(Year)", xlab = "Time(Year)",
        main='RC2')

zz=plot(wavelet3,  plot.cb = TRUE,lwd.coi = 1,col.coi="white",alpha.coi=0.5,col.sig="black",
        lwd.sig = 2, ncol = 768, tol = 0.95, plot.phase=F,ylab="Period(Year)", xlab = "Time(Year)",
        main='RC3')
