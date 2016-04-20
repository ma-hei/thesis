## the kiba dataset

require("openxlsx")
kiba_data = read.xlsx("../data/ci400709d_si_002.xlsx",5,startRow = 1)

## the first column holds the CHEMBL IDs of the compounds
compound_IDs = kiba_data[,1]
target_IDs = names(kiba_data[,-1])

dt_mat = kiba_data[,-1]
dt_mat = as.matrix(dt_mat)
class(dt_mat) = "numeric"

## the names of the targets:
names(kiba_data)

length(unique(kiba_data[,1]))
length(names(kiba_data[,-1]))

## make a triplet format
dt_triplet = matrix(0, nrow = length(which(!is.na(dt_mat))), ncol = 3)
dt_triplet[,1] = which(!is.na(dt_mat), arr.ind=T)[,1]
dt_triplet[,2] = which(!is.na(dt_mat), arr.ind=T)[,2]
dt_triplet[,3] = dt_mat[which(!is.na(dt_mat))]

threshold=8
flag = TRUE
while (flag)
{
  # row stage
  tb = table(dt_triplet[,1])
  ind = which(tb<=threshold)
  cat('removing ',length(ind),' drugs with less than ',threshold,' entries\n')
  # nms = as.numeric(names(tb)[ind])
  nms = names(tb)[ind]
  ind = match(dt_triplet[,1],nms,0)
  ind = which(ind==0)
  dt_triplet = dt_triplet[ind,]
  # col stage
  tb = table(dt_triplet[,2])
  if (min(tb)>threshold) {
    flag = FALSE
  } else {
    ind = which(tb<=threshold)
    cat('removing ',length(ind),' targets with less than ',threshold,' entries\n')
    # nms = as.numeric(names(tb)[ind])
    nms = names(tb)[ind]
    ind = match(dt_triplet[,2],nms,0)
    ind = which(ind==0)
    dt_triplet = dt_triplet[ind,]
  }
  # check row
  tb = table(dt_triplet[,1])
  if (min(tb)>threshold)
    flag = FALSE
}

nrow(dt_triplet)/(length(unique(dt_triplet[,2]))*length(unique(dt_triplet[,1])))

dt_compound_IDs = compound_IDs[unique(dt_triplet[,1])] 
dt_target_IDs = target_IDs[unique(dt_triplet[,2])]

################################
## match CHEMBL IDs to pubchem CIDs
require('httr')
df = data.frame(chembl_id=dt_compound_IDs, cid = 0)
for (i in 1:length(dt_compound_IDs)){
  temp = paste("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",dt_compound_IDs[i],sep ="")
  temp = paste(temp, "/cids/TXT",sep="")
  cat(temp,'\n')
  ans = GET(temp)
  cat(status_code(ans),'\n')
  if (status_code(ans)==200){
    cat(content(ans, encoding = "UTF-8"),'\n')
    df[i,2] = strsplit(gsub("[\n]", " ", content(ans, encoding = "UTF-8"))," ")[[1]][1]
  } else {
    cat('CID not found','\n')
    df[i,2] = "undefined"
  }
}
write(as.matrix(df), file = "chembl_cid_dict.rda")
which(df[,2]=="undefined")
length(which(df[,2]=="undefined"))
undefined_ids = df[which(df[,2]=="undefined"),1]
temp = compound_IDs[dt_triplet[,1]] 
inds = which(temp %in% undefined_ids)
dt_triplet = dt_triplet[-inds,]
##

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))

################
## get the compound similarity matrix
write.table(df[which(df[,2]!="undefined"),2], file="compound_cids.txt", sep="\n", col.names = F, row.names = F, quote=F)
compound_simis = read.csv("..//data//tanimoto_cluster_reqid_3208890451625425647.csv")
compound_ids_simis = compound_simis[,1]
compound_simis = compound_simis[,c(-1,-3267)]
dim(compound_simis)

cids = df[match(compound_IDs[unique(dt_triplet[,1])], df[,1]),2]
which(cids%in%compound_ids_simis)
temp = match(cids, compound_ids_simis)

simi_mat = as.matrix(compound_simis[temp, temp])
## test
temp = compound_IDs[dt_triplet[,1]] 
compound_a = 1
compound_b = 3269
simi_mat[compound_a ,compound_b]
df[match(compound_IDs[unique(dt_triplet[,1])][compound_a], df[,1]),2]
df[match(compound_IDs[unique(dt_triplet[,1])][compound_b], df[,1]),2]
##

## get the target similarity matrix
dt_target_IDs = target_IDs[unique(dt_triplet[,2])]
paste(dt_target_IDs, collapse=' ')
##

dt_compounds = compound_IDs[unique(dt_triplet[,1])]
dt_triplet[,1] = match(dt_triplet[,1], unique(dt_triplet[,1]))
dt_triplet[,2] = match(dt_triplet[,2], unique(dt_triplet[,2]))

save(dt_triplet, simi_mat, file="kiba_data.rda")

counter=0
for (i in 1:n_drugs){
  if (length(which(simi_mat[i,]>0.8))>1){
    counter = counter+1
  }
}
counter

################
## the IC50, Ki mix data set

load('../data/mix_dataset.Rda')

n_drugs = length(unique(mix_dataset[,1]))
n_targets = length(unique(mix_dataset[,2]))


