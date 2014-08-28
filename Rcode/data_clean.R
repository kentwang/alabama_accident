##########################################################
#  This script is to combine and clean raw csv files into
#  a .rdata file.
#
#  Author: Ketong Wang
##########################################################

a <- read.csv("data/accident.csv")
classInfo <- read.csv("data/classInfo.csv")
coordinate <- read.csv("data/coordinate.csv")
