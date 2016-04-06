require("openxlsx")

kiba_data = read.xlsx("../data/ci400709d_si_002.xlsx",5,startRow = 1)

## the first column holds the CHEMBL IDs of the compounds
kiba_data[,1]

## the names of the targets:
names(kiba_data)
