data.all <- read.csv('data/RawData.csv')

data.test <- subset(data.all, experiment==4 | experiment==3)
data.test$experiment <- NULL

# separate target and distractor
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
data.test$stimuli <- trim(data.test$stimuli)
get.target <- function (x) gsub("^t| +.+$", "", x)
get.distractor <- function (x) gsub("^t.+ d", "", x)
data.test$target <- get.target(data.test$stimuli)
data.test$distractor <- get.distractor(data.test$stimuli)
data.test$stimuli <- NULL

# renumber participants
data.test$participant <- as.numeric(factor(data.test$participant))

# remove whitespace
data.test$trial_type <- trim(data.test$trial_type)

write.csv(data.test, file="data/data-merged.csv", row.names = F)
