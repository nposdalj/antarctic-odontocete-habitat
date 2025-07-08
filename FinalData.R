# Script to add all environmental data to daily presence dataframe for final modeling
library(tidyverse)

# ---------------- Step 0: Create base final dataframe-------------
dailyDetection <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Data/dailyDetections.csv")
allData <- dailyDetection
allData$date <- allData$Day
allData <- allData %>% subset(select = -Day)

# --------------- Step 1: Format/Add Antarctic Oscillation Index ----------------
AAO <- read.csv("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Environmental Data/Daily_AAO.csv")

# Make date column
# Add zero in front of single digit dates
addZero <- function(number) {
  if(number == '1') { return('01')
  } else if(number == '2') { return('02')
  } else if(number == '3') { return('03')
  } else if(number == '4') { return('04')
  } else if(number == '5') { return('05')
  } else if(number == '6') { return('06')
  } else if(number == '7') { return('07')
  } else if(number == '8') { return('08')
  } else if(number == '9') { return('09')}
}
# Parse through relevant years to add zero
AAO <- filter(AAO, year >= 2014 & year <= 2016)
AAO$month <- as.character(AAO$month)
AAO$day <- as.character(AAO$day)
for(x in 1:1096) {
  if(AAO$month[x] == '1' | AAO$month[x] == '2' | AAO$month[x] == '3' | AAO$month[x] == '4' | AAO$month[x] == '5' | 
     AAO$month[x] == '6' | AAO$month[x] == '7' | AAO$month[x] == '8' | AAO$month[x] == '9') {
    AAO$month[x] <- addZero(AAO$month[x])
  }
  if(AAO$day[x] == '1' | AAO$day[x] == '2' | AAO$day[x] == '3' | AAO$day[x] == '4' | AAO$day[x] == '5' | 
     AAO$day[x] == '6' | AAO$day[x] == '7' | AAO$day[x] == '8' | AAO$day[x] == '9') {
    AAO$day[x] <- addZero(AAO$day[x])
  }
}
# Create date column with year, month, day
AAO$date <- paste(AAO$year, AAO$month, AAO$day, sep = "")
AAO$date <- as.Date(AAO$date, "%Y%m%d")
# Keep date column & AAO index, only for relevant date ranges
AAO <- AAO %>% subset(select = c(aao_index_cdas, date)) %>%
  filter((date >= '2014-03-05' & date <= '2014-07-17') | 
                        (date >= '2015-02-10' & date <= '2016-01-29') | 
                        (date >= '2016-02-04' & date <= '2016-12-02'))
# Add AAO to final dataset
allData <- merge(allData, AAO, by=intersect(names(allData), names(AAO)))