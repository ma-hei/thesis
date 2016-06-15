args = commandArgs(trailingOnly = T)
file_a = as.numeric(as.numeric(args[1]))
file_b = as.numeric(as.numeric(args[2]))

load(file_a)
load(file_b)

labels = 