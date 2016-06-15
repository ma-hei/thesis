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

save(dt_mat, dt_triplet, compound_IDs, target_IDs, file="kiba_data_raw.rda")

load('kiba_data_raw.rda')

dim(dt_mat)
dim(dt_triplet)

threshold=10
flag = TRUE
while (flag)
{
  # row stage
  tb = table(dt_triplet[,1])
  ind = which(tb<=threshold | tb>200)
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

length(unique(dt_triplet[,1]))
length(unique(dt_triplet[,2]))
nrow(dt_triplet)/(length(unique(dt_triplet[,2]))*length(unique(dt_triplet[,1])))

dt_compound_IDs = compound_IDs[unique(dt_triplet[,1])]


################################
## match CHEMBL IDs to pubchem CIDs to get compound similarities.. 
require('httr')
df = data.frame(chembl_id=dt_compound_IDs, cid = NA)
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
save(df, file = "chembl_cid_dict.rda")
load('chembl_cid_dict.rda')
which(df[,2]=="undefined")
length(which(df[,2]=="undefined"))

undefined_ids = df[which(df[,2]=="undefined"),1]
temp = compound_IDs[dt_triplet[,1]] 
inds = which(temp %in% undefined_ids)
dt_triplet = dt_triplet[-inds,]
##

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))

## At this point we have to compute the drug similarities!!!

## get the target similarity matrix
dt_target_IDs = target_IDs[unique(dt_triplet[,2])]
paste(dt_target_IDs, collapse=' ')
##

dt_triplet[,1] = match(dt_triplet[,1], unique(dt_triplet[,1]))
dt_triplet[,2] = match(dt_triplet[,2], unique(dt_triplet[,2]))

save(dt_triplet, file="kiba_2_triplet.rda")

################

save(sim_mat, file = 'kiba_2_target_sim.rda')

load('kiba_2_triplet.rda')
load('kiba_drug_sim.rda')
load('kiba_target_sim.rda')
write.table(sim_mat, file='kiba_2_target_sim.txt', sep = ' ', col.names = F, row.names = F)
write.table(drug_sim_mat, file='kiba_2_drug_sim.txt', sep=' ', col.names = F, row.names = F)

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))
dt_mat = matrix(NaN, nrow = n_drugs, ncol = n_targets)
dt_mat[dt_triplet[,c(1,2)]] = dt_triplet[,3]
write.table(dt_mat, file='kiba_2_sparse_dt_mat.txt', sep=' ',col.names = F, row.names = F, na = "NaN")

bin_dt_mat = matrix(0, nrow = n_drugs, ncol = n_targets)
bin_dt_mat[dt_triplet[,c(1,2)]] = dt_triplet[,3]
which(bin_dt_mat<=3)
bin_dt_mat[which(bin_dt_mat>3)] = 0
bin_dt_mat[which(bin_dt_mat!=0)] = 1

write.table(bin_dt_mat, file='kiba_2_bin_dt_mat.txt', sep=' ',col.names = F, row.names = F, na = "NaN")



