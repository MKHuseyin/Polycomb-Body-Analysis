# Analysis script for collecting 3D Objects Counter outputs and processing them to extract values and plot

library(readtext)
library(ggplot2)
library(ggsignif)
library(plyr)
library(reshape2)

# before running, add to each stack a measurement made manually in ImageJ of a region of each z-stack with only background signal for normalisation. Put these inside the z-stack folders and name DarkPix.txt

path <- "C:/Users/newto/OneDrive/Documents/Work/Data/DPhil/Microscopy/SD PcG Body analysis/T6/Grouped Files"       # path to directory containing folder containing individual stack folders from ImageJ macros

# setting up

bad.cell <- 0 # to record any cells which did not segment correctly
condition.names <- list.dirs(path, full.names = TRUE, recursive = FALSE) # if directory contains multiple folders containing stack folders, uses their names as conditions (eg UNT and TAM)
conditions <- substr(condition.names, nchar(path) + 2, 200)
allBodies <- data.frame()
allBodies.corrected <- data.frame()
allNuclei <- data.frame()
allDark <- data.frame()
per.cell <- data.frame()
per.cell.corrected <- data.frame() # blank dfs to put data into

# this is the main chunk of reading in the data, preparing to process

for (l in 1:length(condition.names)) {
  
  bodies <- data.frame() # to record all Polycomb body data for the condition
  nuclei <- data.frame() # same with nuclei
  dark <- data.frame(MeanInt = numeric(1), Stack = numeric(1)) # same with background
  mean.dark <- data.frame(MeanInt = numeric(), Stack = numeric()) # means of background stacks
  cellCount <- 1    # variable to count number of cells used
  bodiesOK <- 1     # binary switch used for filtering out nuclei with insufficient bodies
  stack.number <- 0  # iterating stack numbers to write so that the correct DarkPix file can be subtracted
  folder.names <- list.dirs(condition.names[l], full.names = TRUE, recursive = FALSE) # get a list of stack folders
  
  for (k in 1:length(folder.names)) { # loop through all stacks
    file.names <- dir(folder.names[k], pattern =".txt") # get all the text file names (don't care about the images)
    folder <- folder.names[k]
    stack.number <- stack.number + 1
    
    for (i in 1:length(file.names)) { # loop through all files for the stack
      if (grepl("bodies", file.names[i])) { # identifies txt files containing body data to read in
        file <- read.delim(paste(folder, "/", file.names[i], sep=""), header = TRUE, sep = "", dec = ".") # paste() merges strings for full path 
        if (nrow(file)>0) { # can change this number to be more restrictive
          bodiesOK <- 1 # marks this to include downstream
          file$CellNumber <- cellCount # record which cell out of all counted so far this is
          file$StackNumber <- stack.number # record the stack from which the cell came
          bodies <- rbind(bodies, file) # add to bodies df
          
        } else {
          bodiesOK <- 0 # if did not meet requirements, mark the cell as not usable
          bad.cell <- bad.cell + 1
          print(conditions[l])
        }
      } else if (grepl("nucleus", file.names[i])) {
        if (bodiesOK == 1) { # as long as the cell is usable...
          file <- read.delim(paste(folder, "/", file.names[i], sep=""), header = TRUE, sep = "", dec = ".") # paste() merges strings for full path
          file$CellNumber = cellCount
          file$StackNumber <- stack.number
          cellCount <- cellCount + 1
          nuclei <- rbind(nuclei, file) # add nucleus info to nuclei df
        }
      } else if (grepl("DarkPix", file.names[i])) { # also get the background stack
        file <- read.delim(paste(folder, "/", file.names[i], sep=""), header = TRUE, sep = "", dec = ".") # paste() merges strings for full path
        dark.mean <- mean(file$Mean1) # get the mean background signal for the stack
        dark$MeanInt <- dark.mean
        dark$Stack <- stack.number # record stack number and mean value 
        mean.dark <- rbind(mean.dark, dark) # attach to main bg df
        dark <- data.frame(MeanInt = numeric(1), Stack = numeric(1))
        
      }
    }
  }
  
  condition.dark <- data.frame(MeanIntensity = as.numeric(as.character(mean.dark$MeanInt)), Stack = as.numeric(as.character(mean.dark$Stack)), Condition = conditions[1]) # df of all mean background intensities for the condition
  
  
  condition.bodies <- data.frame(Volume = as.numeric(as.character(bodies$Volume)), MeanIntensity =   as.numeric(as.character(bodies$of)),CellNumber = as.numeric(as.character(bodies$CellNumber)), Stack = as.numeric(as.character(bodies$StackNumber)), Condition = conditions[l]) # make a new df with only the useful columns
  
  for (n in 1:stack.number) {
    condition.bodies$MeanIntensity[condition.bodies$Stack == n] <- condition.bodies$MeanIntensity[condition.bodies$Stack == n] - condition.dark$MeanIntensity[condition.dark$Stack == n]  # subtract mean background value from the same stack
  }
  
  condition.bodies$TotalInt <- condition.bodies$MeanIntensity * condition.bodies$Volume  # calculate total intensity for each body by multiplying volume by mean intensity
  condition.bodies.corrected <- condition.bodies[ which(condition.bodies$Volume < 0.5), ]  # filter for bodies larger than 0.5um^3 (to get rid of oversized objects - likely the result of bad segmentation - down to anything with diameter of 1um)
  
  
  condition.nuclei <- data.frame(Volume = as.numeric(as.character(nuclei$Volume)),MeanIntensity = as.numeric(as.character(nuclei$of)),CellNumber = as.numeric(as.character(nuclei$CellNumber)), Stack = as.numeric(as.character(nuclei$StackNumber)), Condition = conditions[l]) # make a new df with only the useful columns
  
  for (n in 1:stack.number) {
    condition.nuclei$MeanIntensity[condition.nuclei$Stack == n] <- condition.nuclei$MeanIntensity[condition.nuclei$Stack == n] - condition.dark$MeanIntensity[condition.dark$Stack == n]  # subtract mean background value from the same stack
  }
  
  condition.nuclei$TotalInt <- condition.nuclei$MeanIntensity * condition.nuclei$Volume  # calculate total intensity for each nucleus by multiplying volume by mean intensity
  
  nrPerCell <- as.data.frame(table(condition.bodies$CellNumber)) # count how many occurrences there are of Polycomb bodies for each cell in new df - useful for a distribution
  nrPerCell$Condition <- conditions[l] # record condition
  
  nrPerCell.corrected <- as.data.frame(table(condition.bodies.corrected$CellNumber)) # count how many occurrences there are of Polycomb bodies for each cell in new df - useful for a distribution
  nrPerCell.corrected$Condition <- conditions[l]
  
  allBodies.corrected <- rbind(allBodies.corrected, condition.bodies.corrected)
  allBodies <- rbind(allBodies, condition.bodies)
  allNuclei <- rbind(allNuclei, condition.nuclei)
  per.cell <- rbind(per.cell, nrPerCell)
  per.cell.corrected <- rbind(per.cell.corrected, nrPerCell.corrected) # update overall dfs for all conditions with data from this one
  
  assign(paste(conditions[l], ".meanvolumes", sep = ""), mean(condition.bodies.corrected$Volume))
  assign(paste(conditions[l], ".meanintensities", sep = ""), mean(condition.bodies.corrected$MeanIntensity))
  assign(paste(conditions[l], ".meanpercell", sep = ""), mean(nrPerCell.corrected$Freq))
  
  assign(paste(conditions[l], ".volumes", sep = ""), condition.bodies.corrected$Volume)
  assign(paste(conditions[l], ".intensities", sep = ""), condition.bodies.corrected$MeanIntensity)
  assign(paste(conditions[l], ".percell", sep = ""), nrPerCell.corrected$Freq)
  assign(paste(conditions[l], ".volumes", sep = ""), condition.bodies.corrected$Volume) # make some quick reference values to check data with
  
}

########### FRACTIONS AND CELL MEANS ##############

# make dfs to fill
allBodies.corrected.fraction <- data.frame(IntFraction = vector(mode = "numeric", length = length(allNuclei$CellNumber)), VolFraction = vector(mode = "numeric", length = length(allNuclei$CellNumber)), Enrichment = vector(mode = "numeric", length = length(allNuclei$CellNumber)), Condition = allNuclei$Condition, CellNumber = allNuclei$CellNumber)
allBodies.corrected.cellmeans <- data.frame(NucleoMeanInt = vector(mode = "numeric", length = length(allNuclei$CellNumber)), MeanBodyMeanInt = vector(mode = "numeric", length = length(allNuclei$CellNumber)), Enrichment = vector(mode = "numeric", length = length(allNuclei$CellNumber)), MeanBodyVol = vector(mode = "numeric", length = length(allNuclei$CellNumber)), MeanBodyTotalInt = vector(mode = "numeric", length = length(allNuclei$CellNumber)), Condition = allNuclei$Condition, CellNumber = allNuclei$CellNumber)

cell.count <- 1

for (j in 1:length(condition.names)) {
  all.condition.bodies.corrected <- allBodies.corrected[ which(allBodies.corrected$Condition == conditions[j]), ] # filter for only ones with same condition
  all.condition.nucleus <- allNuclei[ which(allNuclei$Condition == conditions[j]), ] # same for nuclei
  for (m in 1:max(all.condition.nucleus$CellNumber)) {
    cell.bodies.corrected <- all.condition.bodies.corrected[ which(all.condition.bodies.corrected$CellNumber == m), ]  # take all bodies with matching cell number
    cell.nucleus <- all.condition.nucleus[ which(all.condition.nucleus$CellNumber == m), ] # same for nucleus
    total.intensity.body.corrected <- sum(cell.bodies.corrected$TotalInt) # total intensity from bodies
    cell.body.fraction.int.corrected <- total.intensity.body.corrected / cell.nucleus$TotalInt # calculate the fraction of total nuclear signal arising from bodies
    allBodies.corrected.fraction$IntFraction[cell.count] <- cell.body.fraction.int.corrected # record the fraction in df
    
    total.vol.body.corrected <- sum(cell.bodies.corrected$Volume)
    cell.body.fraction.vol.corrected <- total.vol.body.corrected / cell.nucleus$Volume
    allBodies.corrected.fraction$VolFraction[cell.count] <- cell.body.fraction.vol.corrected # calculate and record fraction of total nuclear volume
    
    mean.meanintensity.body.corrected <- mean(cell.bodies.corrected$MeanIntensity) # calculate means of different metrics and record
    mean.totalintensity.body.corrected <- mean(cell.bodies.corrected$TotalInt)
    mean.volume.body.corrected <- mean(cell.bodies.corrected$Volume)
    allBodies.corrected.cellmeans$MeanBodyMeanInt[cell.count] <- mean.meanintensity.body.corrected
    allBodies.corrected.cellmeans$MeanBodyVol[cell.count] <- mean.volume.body.corrected
    allBodies.corrected.cellmeans$MeanBodyTotalInt[cell.count] <- mean.totalintensity.body.corrected
    
    nucleo.int.corrected <- cell.nucleus$TotalInt - total.intensity.body.corrected # calculate the intensity of each nucleus excluding Polycomb bodies
    nucleo.vol.corrected <- cell.nucleus$Volume - total.vol.body.corrected # repeat for volume
    nucleo.mean.int.corrected <- nucleo.int.corrected / nucleo.vol.corrected # mean intensity for non-Polycomb body volume
    allBodies.corrected.cellmeans$NucleoMeanInt[cell.count] <- nucleo.mean.int.corrected

    cell.count <- cell.count + 1
  }
}

# more general calculations of useful values from above (NB these are for all bodies, will average across conditions there are >1)

mean.volume.body.corrected <- mean(allBodies.corrected$Volume)

mean.nuclear.int.fract.corr <- mean(allBodies.corrected.fraction$IntFraction)
mean.nuclear.vol.fract.corr <- mean(allBodies.corrected.fraction$VolFraction)

################ PREPARING THE DFS #############

# organising DFs for downstream analysis and plotting when >1 condition. Change bracketed values as appropriate. Example for UNT and dTAG conditions

allBodies.corrected$Condition <- factor(allBodies.corrected$Condition, levels = c("UNT", "dTAG"))
per.cell.corrected$Condition <- factor(per.cell.corrected$Condition, levels = c("UNT", "dTAG"))
allNuclei$Condition <- factor(allNuclei$Condition, levels = c("UNT", "dTAG"))
allBodies.corrected.fraction$Condition <- factor(allBodies.corrected.fraction$Condition, levels = c("UNT", "dTAG"))
allBodies.corrected.cellmeans$Condition <- factor(allBodies.corrected.cellmeans$Condition, levels = c("UNT", "dTAG"))

# separate conditions into new dfs. Add more lines as required.

UNT.allBodies.fraction <- allBodies.corrected.fraction[which(allBodies.corrected.fraction$Condition == "UNT"), ]
dTAG.allBodies.fraction <- allBodies.corrected.fraction[which(allBodies.corrected.fraction$Condition == "dTAG"), ]


########### VOLUME QUARTILES ###########
# (for multiple conditions)

# split all bodies into individual conditions
dTAG.allBodies.corrected = allBodies.corrected[ which(allBodies.corrected$Condition == "dTAG"), ]
UNT.allBodies.corrected = allBodies.corrected[ which(allBodies.corrected$Condition == "UNT"), ]

# add column for quartile variable
UNT.allBodies.corrected$vol_quartiles = NA
dTAG.allBodies.corrected$vol_quartiles = NA

# define quartiles based on UNT
UNT.allBodies.corrected$vol_quartiles = cut(UNT.allBodies.corrected$Volume, quantile(UNT.allBodies.corrected$Volume, probs = seq(from = 0,to = 1, length.out = 5), type = 5, names = FALSE), labels = FALSE, include.lowest=TRUE)

# make quartiles into factors
UNT.allBodies.corrected$vol_quartiles <- as.factor(UNT.allBodies.corrected$vol_quartiles)

UNT.allBodies.corrected$vol_quartiles <- mapvalues(UNT.allBodies.corrected$vol_quartiles, from = c('1', '2', '3', '4'), to = c('Q1', 'Q2', 'Q3', 'Q4'))

# determine quartile bounds in UNT
Q4.min <- min(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$vol_quartiles == "Q4"), "Volume"])
Q3.max <- max(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$vol_quartiles == "Q3"), "Volume"])
Q3.min <- min(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$vol_quartiles == "Q3"), "Volume"])
Q2.max <- max(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$vol_quartiles == "Q2"), "Volume"])
Q2.min <- min(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$vol_quartiles == "Q2"), "Volume"])
Q1.max <- max(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$vol_quartiles == "Q1"), "Volume"])

# assign treated volumes to quartiles based on UNT boundaries
for (m in 1:length(dTAG.allBodies.corrected$Volume)) {
  if (dTAG.allBodies.corrected$Volume[m] >= Q4.min) {
    dTAG.allBodies.corrected$vol_quartiles[m] <- "Q4"
  } else if (dTAG.allBodies.corrected$Volume[m] >= Q3.min && dTAG.allBodies.corrected$Volume[m] <= Q3.max) {
    dTAG.allBodies.corrected$vol_quartiles[m] <- "Q3"
  } else if (dTAG.allBodies.corrected$Volume[m] >= Q2.min && dTAG.allBodies.corrected$Volume[m] <= Q2.max) {
    dTAG.allBodies.corrected$vol_quartiles[m] <- "Q2"
  } else if (dTAG.allBodies.corrected$Volume[m] <= Q1.max) {
    dTAG.allBodies.corrected$vol_quartiles[m] <- "Q1"
  }
}

dTAG.allBodies.corrected$vol_quartiles <- as.factor(dTAG.allBodies.corrected$vol_quartiles)

# make a new column for counts per quartile
per.cell.corrected$Q1ct <- NA
per.cell.corrected$Q2ct <- NA
per.cell.corrected$Q3ct <- NA
per.cell.corrected$Q4ct <- NA

# count how many are in each quartile in each cell
per.cell.corrected$Var1 <- as.numeric(per.cell.corrected$Var1)
for (n in 1:max(per.cell.corrected[ which(per.cell.corrected$Condition == "UNT"), "Var1"])) {
  per.cell.corrected[ which(per.cell.corrected$Condition == "UNT" & per.cell.corrected$Var1 == n), "Q1ct"] <- length(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$CellNumber == n & UNT.allBodies.corrected$vol_quartiles == "Q1"), "CellNumber"])
  per.cell.corrected[ which(per.cell.corrected$Condition == "UNT" & per.cell.corrected$Var1 == n), "Q2ct"] <- length(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$CellNumber == n & UNT.allBodies.corrected$vol_quartiles == "Q2"), "CellNumber"])
  per.cell.corrected[ which(per.cell.corrected$Condition == "UNT" & per.cell.corrected$Var1 == n), "Q3ct"] <- length(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$CellNumber == n & UNT.allBodies.corrected$vol_quartiles == "Q3"), "CellNumber"])
  per.cell.corrected[ which(per.cell.corrected$Condition == "UNT" & per.cell.corrected$Var1 == n), "Q4ct"] <- length(UNT.allBodies.corrected[ which(UNT.allBodies.corrected$CellNumber == n & UNT.allBodies.corrected$vol_quartiles == "Q4"), "CellNumber"])
}
for (n in 1:max(per.cell.corrected[ which(per.cell.corrected$Condition == "dTAG"), "Var1"])) {
  per.cell.corrected[ which(per.cell.corrected$Condition == "dTAG" & per.cell.corrected$Var1 == n), "Q1ct"] <- length(dTAG.allBodies.corrected[ which(dTAG.allBodies.corrected$CellNumber == n & dTAG.allBodies.corrected$vol_quartiles == "Q1"), "CellNumber"])
  per.cell.corrected[ which(per.cell.corrected$Condition == "dTAG" & per.cell.corrected$Var1 == n), "Q2ct"] <- length(dTAG.allBodies.corrected[ which(dTAG.allBodies.corrected$CellNumber == n & dTAG.allBodies.corrected$vol_quartiles == "Q2"), "CellNumber"])
  per.cell.corrected[ which(per.cell.corrected$Condition == "dTAG" & per.cell.corrected$Var1 == n), "Q3ct"] <- length(dTAG.allBodies.corrected[ which(dTAG.allBodies.corrected$CellNumber == n & dTAG.allBodies.corrected$vol_quartiles == "Q3"), "CellNumber"])
  per.cell.corrected[ which(per.cell.corrected$Condition == "dTAG" & per.cell.corrected$Var1 == n), "Q4ct"] <- length(dTAG.allBodies.corrected[ which(dTAG.allBodies.corrected$CellNumber == n & dTAG.allBodies.corrected$vol_quartiles == "Q4"), "CellNumber"])
}

# melt into plot-able dfs
plot.per.cell.corrected <- melt(per.cell.corrected, id.vars = "Condition", measure.vars = c("Q1ct", "Q2ct", "Q3ct", "Q4ct"))

plot.per.cell.corrected$variable <- mapvalues(plot.per.cell.corrected$variable, from = c("Q1ct", "Q2ct", "Q3ct", "Q4ct"), to = c("Q1", "Q2", "Q3", "Q4"))

########### ENRICHMENT CALCS ###########

plot.allBodies.cellmeans <- melt(allBodies.corrected.cellmeans, id.vars = "Condition", measure.vars = c("NucleoMeanInt", "MeanBodyMeanInt"))
plot.allBodies.cellmeans$variable <- mapvalues(plot.allBodies.cellmeans$variable, from = c("NucleoMeanInt", "MeanBodyMeanInt"), to = c("Nucleoplasm", "Polycomb Body"))
plot.allBodies.cellmeans$variable <- factor(plot.allBodies.cellmeans$variable, levels = c("Polycomb Body", "Nucleoplasm")) # order so bodies first
plot.allBodies.cellmeans$Condition <- factor(plot.allBodies.cellmeans$Condition, levels = c("UNT", "dTAG"))

UNT.median.nucleo <- median(UNT.cellmeans$NucleoMeanInt) # measure median for non-body for each condition (so that median line is aligned with 1)
dTAG.median.nucleo <- median(dTAG.cellmeans$NucleoMeanInt)
norm.plot.allBodies.cellmeans <- plot.allBodies.cellmeans
norm.plot.allBodies.cellmeans$value <- ifelse(norm.plot.allBodies.cellmeans$Condition == "UNT", norm.plot.allBodies.cellmeans$value / UNT.median.nucleo, norm.plot.allBodies.cellmeans$value) # normalise each condition's mean Polycomb body signal density with the corresponding median
norm.plot.allBodies.cellmeans$value <- ifelse(norm.plot.allBodies.cellmeans$Condition == "dTAG", norm.plot.allBodies.cellmeans$value / dTAG.median.nucleo, norm.plot.allBodies.cellmeans$value)

######################################### PLOTS #####################################################
# there are lots of options to plot here. These are the most relevant.

############# Nuclear volumes (to check for any changes in gross nuclear morphology)
p <- ggplot(data = allNuclei, aes(x = Condition, y = Volume, fill = Condition)) +  
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "red"), breaks=c("UNT", "dTAG")) +
  geom_signif(comparisons = list(c(1, 2)),
              step_increase = 0.1,
              test = "wilcox.test") +
  expand_limits(y = 0)

############### Nucleus total intensity (check for changes in protein expression)
p <- ggplot(data = allNuclei, aes(x = Condition, y = TotalInt, fill = Condition)) +  # nuclei total intensity boxplot
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "red"), breaks=c("UNT", "dTAG")) +
  geom_signif(comparisons = list(c(1, 2)),
              step_increase = 0.1,
              test = "t.test") +
  expand_limits(y = 0)

##############  Body volume
p <- ggplot(data = allBodies.corrected, aes(x = Condition, y = Volume, fill = Condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "red"), breaks=c("UNT", "dTAG")) +
  geom_signif(comparisons = list(c(1, 2)),
              step_increase = 0.1,
              test = "wilcox.test") +
  expand_limits(y = 0)

##############  Body mInt (signal density from bodies)
p <- ggplot(data = allBodies.corrected, aes(x = Condition, y = MeanIntensity, fill = Condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "red"), breaks=c("UNT", "dTAG")) +
  geom_signif(comparisons = list(c(1, 2)),
              step_increase = 0.1,
              test = "wilcox.test") +
  expand_limits(y = 0)

##############  per cell (numbers of bodies)
p <- ggplot(data = per.cell.corrected, aes(x = Condition, y = Freq, fill = Condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "red"), breaks=c("UNT", "dTAG")) +
  geom_signif(comparisons = list(c(1, 2)),
              step_increase = 0.1,
              test = "t.test",
              y_position = 300) +
  ylim(0, 300)

##############  Body volume quartiles
# plot all together (no stats)
p <- ggplot(data = plot.per.cell.corrected, aes(x = variable, y = value, fill = Condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "red"), breaks=c("UNT", "dTAG")) +
  ylab("Bodies per cell") +
  expand_limits(y = 0)
print(p)
dev.off()

# plot in paired panels - allows stats
p <- ggplot(data = plot.per.cell.corrected, aes(x = Condition, y = value, fill = Condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "red"), breaks=c("UNT", "dTAG")) +
  ylab("Bodies per cell") +
  expand_limits(y = 0) +
  facet_wrap(vars(variable), scales = "free_x") +
  geom_signif(comparisons = list(c(1, 2)),
              step_increase = 0.1,
              test = "wilcox.test")

############## nucleo mean int vs body mean int (enrichment of body relative to non-body)
# plot together (no stats)
p <- ggplot(data = norm.plot.allBodies.cellmeans, aes(x = variable, y = value, fill = Condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "green", "purple"), breaks=c("UNT", "dTAG", "dTAGTAM")) +
  ylab("Mean Intensity")

# plot individual panels - stats
p <- ggplot(data = norm.plot.allBodies.cellmeans, aes(x = Condition, y = value, fill = Condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "green", "purple"), breaks=c("UNT", "dTAG", "dTAGTAM")) +
  ylab("Bodies per cell") +
  expand_limits(y = 0) +
  facet_wrap(vars(variable), scales = "free_x") +
  geom_signif(comparisons = list(c(1, 2)),
              step_increase = 0.1,
              test = "wilcox.test")