
# go where the data is
setwd('~/Documents/Research/FG2/analysis/')

library("robustbase")
library("MASS")
library("car")
library("reshape2")
library("ggplot2")
library('lme4')
library('lmerTest')

#########################
# read and clean the data
#########################
dat <- read.csv('~/Documents/Research/FG2/analysis/Cond_ROI_05152020.csv')
names(dat)[1] <- 'ID'
names(dat) <- gsub('Hip', 'HIP', names(dat))
names(dat) <- gsub('Amy', 'AMY', names(dat))

# no specifier for lh+rh ROIs 
names(dat) <- gsub('Comb', '', names(dat), ignore.case = TRUE)

stria <- read.csv('~/Documents/Research/FG2/analysis/Cond_striatum_hemi_combined.csv')
names(stria)[1] <- 'ID'

dat <- merge(x=dat, y=stria, by = 'ID')
names(dat) <- gsub('executive', 'EXE', names(dat))
names(dat) <- gsub('limbic', 'LIM', names(dat))
names(dat) <- gsub('insula', 'INS', names(dat), ignore.case = TRUE)
names(dat) <- gsub('AngularGyrus', 'ANG', names(dat), ignore.case = TRUE)
# used ROIs 

ROI <- c("EXE","THA", "LIM", "HIP", "AMY", 
         "INS","dmPFC", "PCC","ANG" )

vars <- c('ID', 'Pers.x', 'PDI.x', paste(ROI, '_CSp', sep=''), paste(ROI, '_CSm', sep='')) 
dat <- dat[, vars]

# calculate contrasts 
for (i in 1:length(ROI)){
  dat[[paste(ROI[i], "contrast", sep='')]] <- dat[[paste(ROI[i], "_CSp", sep='')]] - dat[[paste(ROI[i], "_CSm", sep='')]]
}

#########################
# Plot main effects 
#########################

# for each ROI, plot CS+ vs CS-
pltdat <- dat[, c('ID', paste(ROI, "_CSp", sep=''), paste(ROI, "_CSm", sep='') )]
pltdat <- melt(pltdat, id.vars = "ID" )
pltdat$ROI <- gsub('.{4}$', '', pltdat$variable)
pltdat$cond <- ifelse(grepl('CSp', pltdat$variable, fixed = TRUE), 'CSp', 'CSm' )
pltdat$cond <- factor(pltdat$cond,levels=c('CSp','CSm'),ordered = TRUE)
pltdat$ROI <- factor(pltdat$ROI, levels = unique(pltdat$ROI))

pvals <- vector()
for(r in ROI){
  m <- t.test(dat[[paste(r, '_CSp', sep='')]], dat[[paste(r, '_CSm', sep='')]], paired = TRUE, alternative = "two.sided")
  pvals <- c(pvals, m$p.value)
}

agg <- aggregate(pltdat[, 'value'], by = list(pltdat$ROI), max)
agg <- agg$x
## plot all effects 

plt <- ggplot(pltdat, aes(fill=cond, y=value, x=ROI)) + 
  geom_boxplot(outlier.shape = NA,alpha = 0.5) +
  scale_fill_manual(values=c( "firebrick1",'dodgerblue2')) + 
  geom_point(position=position_jitterdodge(), size = 0.5, aes(color = cond) ) +
  theme(axis.text.x = element_text(angle = 45, size=12, hjust = 1),
        axis.title.x = element_text(size=16, hjust = 0.5 ),
        axis.title.y = element_text(size=16, hjust = 0.5 ),
        axis.text.y = element_text(size=12, hjust = 0.5 )) +
  scale_color_manual(values = c("CSp" = "firebrick1", "CSm" = 'dodgerblue2')) +
  scale_x_discrete(name ="Regions of Interest")
  
for (i in seq(pvals)){
  p <- pvals[i]
  if (p < 0.05 & p >= 0.001) {
    plt <- plt+annotate(geom="text", x=i, y= agg[i]+0.5, label="*",
                    color="black", size=10)
  }
  if (p < 0.001) {
  plt <- plt+annotate(geom="text", x=i, y=agg[i]+0.5, label="**",
                      color="black", size=10)
  }
}

plt

# Check how correlated the effects are 
pairs(dat[,c(paste(ROI, 'contrast', sep =''))])

##################################################
# Fit robust regression between ROI data and PDI 
##################################################

## compute association between PDI and brain responses 
asso1 <- vector()
tscore1 <- vector()
beta1 <- vector()
for (r in ROI){
  r <- paste(r, 'contrast', sep='')
  s<-summary(lmrob(dat[[r]] ~ dat$PDI))
  asso1 <- c(asso1, s$coefficients[2,4])
  tscore1 <- c(tscore1, s$coefficients[2,3])
  beta1 <- c(beta1, s$coefficients[2,1])
}

rmat <- data.frame(t(matrix(nrow=length(beta1), ncol=3, c(round(beta1,3), round(tscore1,2) , round(asso1,3)))))
names(rmat) <- ROI

## test if SQRT transform changes things  
dat$sqrtPDI <- sqrt(dat$PDI)
asso2 <- vector()
for (r in ROI){
  r <- paste(r, 'contrast', sep='')
  s<-summary(lmrob(dat[[r]] ~ dat$sqrtPDI))
  asso2 <- c(asso2, s$coefficients[2,4])
}


##################################################
# Plot all associations  
##################################################

# make an outputdir 
pltsdir <- '~/Documents/Research/FG2/analysis/plts'
if (!dir.exists(pltsdir)) dir.create(pltsdir)

yax <- c(1.25, 1.5, 1.25 , 1.5, 1.25, 1.3, 1.25, 1.5, 1.5)

for (r in seq(ROI)){
  
  y = dat[[paste(ROI[r], 'contrast', sep= '')]]
  
plt <- ggplot(data = dat, aes(x = PDI.x, y = y )) + 
  geom_point(color='dodgerblue1', shape = 21, size=2)+
  stat_smooth(method=function(formula,data,weights=weight) rlm(formula,
                                                               data,
                                                               weights=weight,
                                                               method="MM"),
              fullrange=TRUE)+
  
  labs(title =ROI[r], x = "PDI", y = 'CS+ vs CS- contrast') +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),
        plot.title = element_text(hjust = 0.5, size=16),legend.text = element_text(size=14)) +
  annotate('text' , x = 8, y = yax[r], label = sprintf('Robust regression: \nb=%.3f, t=%.3f, p=%.3f', rmat[[r]][1], rmat[[r]][2], rmat[[r]][3]), size = 3.5, hjust = 0 )  


jpeg(paste(pltsdir, '/', ROI[r], '_PDI_robust.jpg', sep=''), width = 1200, height = 1800, res = 300, units = 'px')

print(plt)
dev.off()
  }
