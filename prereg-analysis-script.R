#
# Analysis for numerical magnitude visual search experiment
#



## LOAD REQUIRED PACKAGES ####
require(runjags)
require(coda)
require(plyr)
source('plotPost.R')



## LOAD THE DATA FILE ####
data.all <- read.csv('data/data-merged.csv')

# subset of all non-practice data
data.experiment <- subset(data.all, trial_type=='experiment')

## data file has the following columns ####
# participant - numeric ID
# global_trial_num - numeric trial number
# block_num - numeric block number
# trial_type - "practice" or "experiment"
# local_trial_num - numeric trial number within block
# targetpresence - 0 if target absent, 1 if target present
# setsize - 2, 3, 4, or 6
# targetlocation - which of the 12 possible locations the target appeared in
# response - 0 for target absent, 1 for target present, -1 for no response before timeout
# correctness - 0 if incorrect, 1 if correct
# responsetime - in seconds
# target - target stimulus
# distractor - distractor stimulus



## REMOVE SUBJECTS WITH LOW ACCURACY ####
# Anyone who is 2SDs below the group average

# Subject accuracy
acc.subject <- ddply(data.experiment, .(participant), function(s){
  return(c(acc=mean(s$correctness)))
})

# find the cutoff that is 2SDs below the mean
sd.acc <- sd(acc.subject$acc)
mean.acc <- mean(acc.subject$acc)
cutoff.acc <- mean.acc - 2*sd.acc

# mark bad subjects
acc.subject$bad <- acc.subject$acc < cutoff.acc

# get list of bad subjects
bad.subjects <- as.character(subset(acc.subject, bad==T)$participant)

# filter out bad subjects
data.experiment.filtered <- subset(data.experiment, !participant%in%bad.subjects)



## SUMMARIZE DATA AT THE SUBJECT LEVEL ####
data.subject <- ddply(subset(data.experiment.filtered, correctness==1), .(participant, targetpresence, setsize, target, distractor), function(s){
 return(c(meanrt = mean(s$responsetime*1000))) 
})

# set of target present trials for analysis
data.subject.targetpresent <- subset(data.subject, targetpresence==1)

## PREPARE DATA FOR JAGS ####
data.jags <- data.subject.targetpresent

# renumber the subjects to avoid gaps
data.jags$participant <- as.numeric(factor(data.subject.targetpresent$participant))

# convert target and distractor into numeric levels
data.jags$target <- as.numeric(data.jags$target)
data.jags$distractor <- as.numeric(data.jags$distractor)

# calculate sd of cell level data for priors
sd.cells <- sd(data.jags$meanrt)

# get vague prior with mode equal to empirical sd
gammaShRaFromModeSD <- function( mode , sd ) {
  if ( mode <=0 ) stop("mode must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

sd.cells.gamma <- gammaShRaFromModeSD( sd.cells, sd.cells*5 )

# calculate subject effect sd

# first summarize mean rt at the subject level
subject.rt <- ddply(subset(data.experiment.filtered, correctness==1 & targetpresence==1), .(participant), function(s){
  return(c(mean.rt = mean(s$responsetime*1000)))
})

# find sd of subject differences
sd.subject <- sd(subject.rt$mean.rt)

# get gamma priors based on sd of subject differences
sd.subject.gamma <- gammaShRaFromModeSD( sd.subject, sd.subject*5)

# prior for intercept
intercept.prior.mean <- mean(subject.rt$mean.rt)
intercept.prior.sd <- sd.subject * 5

# prior for slope
slope.prior.mean <- 0
slope.prior.sd <- 250

# create JAGS data list
jags.data.list <- list(
  rt = data.jags$meanrt,
  target.type = data.jags$target,
  setsize = data.jags$setsize,
  subject = data.jags$participant,
  ndata = nrow(data.jags),
  nsubjects = length(unique(data.jags$participant)),
  ntargettypes = length(unique(data.jags$target)),
  sd.cells.prior.shape = sd.cells.gamma$shape, 
  sd.cells.prior.rate = sd.cells.gamma$rate,
  sd.subject.prior.shape = sd.subject.gamma$shape, 
  sd.subject.prior.rate = sd.subject.gamma$rate,
  intercept.prior.mean = intercept.prior.mean,
  intercept.prior.sd = intercept.prior.sd,
  slope.prior.mean = slope.prior.mean,
  slope.prior.sd = slope.prior.sd
)

## RUN THE JAGS MODEL ####
params <- c('intercept.group', 'slope.group', 'intercept', 'slope')
results <- run.jags("jags-model-prereg.txt", monitor=params, n.chains=4, data = jags.data.list, adapt=1000,burnin=4000, sample=10000)

plot(results)

#print(results)

codaSamples = as.mcmc.list( results )
mcmcMat = as.matrix(codaSamples)



## CONTRASTS ON POSTERIORS ####

# For reference, the intercept[] and slope[] parameters are estimates of the search function
# at the group level, by different target type. The target types are:
# [1] - Numeric 2
# [2] - Numeric 5
# [3] - Numeric 6
# [4] - Numeric 9
# [5] - Letter b
# [6] - Letter d
# [7] - Flipped 6
# [8] - Flipped 9

# Hypothesis 1: Search slope will be shallower (numerically smaller) 
# with higher magnitude targets compared to lower magnitude targets

slopes.larger.avg = (mcmcMat[,'slope[2]'] + mcmcMat[,'slope[4]'])/2
slopes.smaller.avg = (mcmcMat[,'slope[1]'] + mcmcMat[,'slope[3]'])/2
slopes.magnitude.diff = slopes.larger.avg - slopes.smaller.avg

# The hypothesis predicts that the 95% HDI of slopes.diff will be
# all negative values.
plotPost(slopes.magnitude.diff)

# Hypothesis 2: There will be no differences in slope for letters or
# flipped numbers
slope.letter.diff = mcmcMat[,'slope[5]'] - mcmcMat[,'slope[6]']
plotPost(slope.letter.diff)

slope.flipped.diff = mcmcMat[,'slope[7]'] - mcmcMat[,'slope[8]']
plotPost(slope.flipped.diff)

# Hypothesis 3: The difference in the magnitude should be larger
# than the difference in the letter or flipped conds

# Both of these plots should have HDIs above 0.
plotPost(abs(slopes.magnitude.diff) - abs(slope.letter.diff))
plotPost(abs(slopes.magnitude.diff) - abs(slope.flipped.diff))
