files = list.files()

temp = load(files[1])
temp = as.matrix(temp)
temp = as.vector(temp)
n = length(temp)

predictions = rep(NA, n)

for (i in 1:length(files)){
  temp = load(files[i])
  temp = as.matrix(temp)
  temp = as.vector(temp)
  predictions[which(!is.na(temp))] = temp[which(!is.na(temp))]
}

