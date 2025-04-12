# Alex Searle-Barnes
# 31 January 2024

# Calculate the area under the curve for three peaks in each of the three CT scan images in the manuscript

library(tidyverse)
library(DescTools)
library(reshape)
library(wesanderson)

kv40 <- read.csv(file = "40kV bin 4 0240.csv", header = TRUE)
kv60 <- read.csv(file = "60kV bin 2 0480.csv", header = TRUE)
kv80 <- read.csv(file = "80kV bin 2 0480.csv", header = TRUE)

# Remove first data point
kv40 <- kv40[2:nrow(kv40),]
kv60 <- kv60[2:nrow(kv60),]
kv80 <- kv80[2:nrow(kv80),]

# Identify peaks

# Air
kv40AirInt <- AUC(x = kv40$bin.start, y = kv40$count, from = min(16000), to = max(28500))
kv40AirIntAvg <- kv40AirInt / (28500 - 16000)
kv40AirIntAvg <- kv40AirIntAvg %>% round(digits = 0)
kv40AirMax <- kv40$count[70:122] %>% max()
kv40AirMax2 <- kv40$bin.start[which(kv40$count == kv40AirMax)] 

kv60AirInt <- AUC(x = kv60$bin.start, y = kv40$count, from = min(9000), to = max(15500))
kv60AirIntAvg <- kv60AirInt / (15500 - 9000)
kv60AirIntAvg <- kv60AirIntAvg %>% round(digits = 0)
kv60AirMax <- kv60$count[70:122] %>% max()
kv60AirMax2 <- kv60$bin.start[which(kv60$count == kv60AirMax)] 

kv80AirInt <- AUC(x = kv80$bin.start, y = kv40$count, from = min(9000), to = max(15500))
kv80AirIntAvg <- kv80AirInt / (15500 - 9000)
kv80AirIntAvg <- kv80AirIntAvg %>% round(digits = 0)
kv80AirMax <- kv80$count[70:126] %>% max()
kv80AirMax2 <- kv80$bin.start[which(kv80$count == kv80AirMax)] 


# Foam


kv40FoamInt <- AUC(x = kv40$bin.start, y = kv40$count, from = min(25934), to = max(31162))
kv40FoamIntAvg <- kv40FoamInt / (31162 - 25934)
kv40FoamIntAvg <- kv40FoamIntAvg %>% round(digits = 0)
kv40FoamMax <- kv40 %>% filter(between(bin.start, 25934, 31162)) %>% select(count) %>% max()
kv40FoamMax2 <- kv40$bin.start[which(kv40$count == kv40FoamMax)]

kv60FoamInt <- AUC(x = kv60$bin.start, y = kv40$count, from = min(15112), to = max(18107))
kv60FoamIntAvg <- kv60FoamInt / (18107 - 15112)
kv60FoamIntAvg <- kv60AirIntAvg %>% round(digits = 0)
kv60FoamMax <- kv60 %>% filter(between(bin.start, 15112, 18107)) %>% select(count) %>% max()
kv60FoamMax2 <- kv60$bin.start[which(kv60$count == kv60FoamMax)] 

kv80FoamInt <- AUC(x = kv80$bin.start, y = kv40$count, from = min(14477), to = max(19132))
kv80FoamIntAvg <- kv80FoamInt / (19132 - 14477)
kv80FoamIntAvg <- kv80FoamIntAvg %>% round(digits = 0)
kv80FoamMax <- kv80 %>% filter(between(bin.start, 14477, 19132)) %>% select(count) %>% max()
kv80FoamMax2 <- kv80$bin.start[which(kv80$count == kv80FoamMax)] 



# Calcite
kv40CalcInt <- AUC(x = kv40$bin.start, y = kv40$count, from = min(36795), to = max(58197))
kv40CalcIntAvg <- kv40CalcInt / (58197 - 36795)
kv40CalcIntAvg <- kv40CalcIntAvg %>% round(digits = 0)
kv40CalcMax <- kv40 %>% filter(between(bin.start, 36795, 58197)) %>% select(count) %>% max()
kv40CalcMax2 <- kv40$bin.start[which(kv40$count == kv40CalcMax)] 

kv60CalcInt <- AUC(x = kv60$bin.start, y = kv40$count, from = min(19175), to = max(33341))
kv60CalcIntAvg <- kv60CalcInt / (33341 - 19175)
kv60CalcIntAvg <- kv60CalcIntAvg %>% round(digits = 0)
kv60CalcMax <- kv60 %>% filter(between(bin.start, 19175, 33341)) %>% select(count) %>% max()
kv60CalcMax2 <- kv60$bin.start[which(kv60$count == kv60CalcMax)] 

kv80CalcInt <- AUC(x = kv80$bin.start, y = kv40$count, from = min(20871), to = max(31675))
kv80CalcIntAvg <- kv80CalcInt / (31675 - 20871)
kv80CalcIntAvg <- kv80CalcIntAvg %>% round(digits = 0)
kv80CalcMax <- kv80 %>% filter(between(bin.start, 20871, 31675)) %>% select(count) %>% max()
kv80CalcMax2 <- kv80$bin.start[which(kv80$count == kv80CalcMax)] 

#Calculate distances between peaks
kv40distFoamCalc <- kv40CalcMax2 - kv40FoamMax2
kv60distFoamCalc <- kv60CalcMax2 - kv60FoamMax2
kv80distFoamCalc <- kv80CalcMax2 - kv80FoamMax2

# Calculate calcite to foam count ratio
kv40CalcFoamRatio <- kv40CalcMax / kv40FoamMax
kv60CalcFoamRatio <- kv60CalcMax / kv60FoamMax
kv80CalcFoamRatio <- kv80CalcMax / kv80FoamMax




# Output useful table
numbersOut <- data.frame(Figure = c("A", "B", "C"),
                         kV = c("40", "60", "80"),
                         airPeakCount = c(kv40AirMax, kv60AirMax, kv80AirMax),
                         airPeakBright = c(kv40AirMax2, kv60AirMax2, kv80AirMax2),
                         airIntMean = c(kv40AirIntAvg, kv60AirIntAvg, kv80AirIntAvg),
                         foamPeakCount = c(kv40FoamMax, kv60FoamMax, kv80FoamMax),
                         foamPeakBright = c(kv40FoamMax2, kv60FoamMax2, kv80FoamMax2),
                         foamIntMean = c(kv40FoamIntAvg, kv60FoamIntAvg, kv80FoamIntAvg),
                         calcitePeakCount = c(kv40CalcMax, kv60CalcMax, kv80CalcMax),
                         calcitePeakBright = c(kv40CalcMax2, kv60CalcMax2, kv80CalcMax2),
                         calciteIntMean = c(kv40CalcIntAvg, kv60CalcIntAvg, kv80CalcIntAvg),
                         distFoamCalc = c(kv40distFoamCalc, kv60distFoamCalc, kv80distFoamCalc),
                         CalcFoamRatio = c(kv40CalcFoamRatio, kv60CalcFoamRatio, kv80CalcFoamRatio)
                         )

write.csv(numbersOut, file = "Histogram peak calculations.csv", row.names = FALSE)

# Make some graphs

kv40plot <- ggplot(kv40, aes(x = index, y = count))+
  geom_point()+
  geom_line()+
  labs(title = "40 kV",
       x = "Index",
       y = "Count")

kv60plot <- ggplot(kv60, aes(x = index, y = count))+
  geom_point()+
  geom_line()+
  labs(title = "60 kV",
       x = "Index",
       y = "Count")

kv80plot <- ggplot(kv80, aes(x = index, y = count))+
  geom_point()+
  geom_line()+
  labs(title = "80 kV",
       x = "Index",
       y = "Count")

# combine data
allData <- data.frame(index = kv40$index,
                      kv40bin = kv40$bin.start,
                      kv40count = kv40$count,
                      kv60bin = kv60$bin.start,
                      kv60count = kv60$count,
                      kv80bin = kv80$bin.start,
                      kv80count = kv80$count)


# Melt to make longer by index
allData2 <- allData %>% select(index, kv40count, kv60count, kv80count) %>% melt(., id = "index")


# Make plot
allPlots <- ggplot(allData)+
  geom_point(aes(x = kv40bin, y = kv40count, colour = "40 kV bin 4"))+
  geom_line(aes(x = kv40bin, y = kv40count, colour = "40 kV bin 4"))+
  geom_point(aes(x = kv60bin, y = kv60count, colour = "60 kV bin 4"))+
  geom_line(aes(x = kv60bin, y = kv60count, colour = "60 kV bin 4"))+
  geom_point(aes(x = kv80bin, y = kv80count, colour = "80 kV bin 4"))+
  geom_line(aes(x = kv80bin, y = kv80count, colour = "80 kV bin 4"))+
  scale_x_continuous(breaks=seq(0, 60000, 10000))+
  scale_colour_manual(labels = c("40 kV bin 4", "60 kV bin2", "80 kV bin 2"), values = wes_palette("Darjeeling1"))+
  labs(title = "X-ray CT brightness values of a foraminifer \nembedded in foam within a plastic straw at increasing source powers",
       subtitle = "",
       x = "Brightness",
       y = "Count",
       color = "Source power")+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave(filename = "allPlots.png", allPlots, width = 5000, height = 5000, dpi = 500, units = "px")

