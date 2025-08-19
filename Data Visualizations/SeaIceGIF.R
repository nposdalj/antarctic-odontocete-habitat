library(tidyverse)
library(magick)

# GIF of sea ice changes (monthly) over three sites
# intervals: 02/15 to 12/16 (KGI and CI data)
#    skipping EI because less ice and data would not be continuous
#    images taken from the first of each month
# images from: https://data.seaice.uni-bremen.de/databrowser/


# adding images to list for GIF
frames <- paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Presentation Plots/Ice Images/", 
                      2:12, "-01-2015.png")
frames <- append(frames, paste0("C:/Users/HARP/Documents/GitHub/antarctic-odontocete-habitat/Presentation Plots/Ice Images/", 
                                          1:12, "-01-2016.png"))

# making GIF
gif <- image_read(frames)
gif <- image_animate(gif, fps=2)
