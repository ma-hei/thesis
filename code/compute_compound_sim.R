load('chembl_cid_dict.rda')
################
## get the compound similarity matrix
which(duplicated(df[which(df[,2]!='undefined'),2]))

#write.table(df[which(df[,2]!="undefined"),2], file="compound_cids.txt", sep="\n", col.names = F, row.names = F, quote=F)
temp = compound_IDs[unique(dt_triplet[,1])] 
write.table(df[match(temp, df[,1]),2], file="compound_cids.txt", sep="\n", col.names = F, row.names = F, quote=F)
compound_sim = read.csv("tanimoto_cluster_reqid_2439928264593734584.csv")
compound_ids_sim = compound_sim[,1]
compound_sim = compound_sim[,c(-1,-2115)]
dim(compound_sim)

cids = df[match(compound_IDs[unique(dt_triplet[,1])], df[,1]),2]
length(which(cids%in%compound_ids_sim))
temp = match(cids, compound_ids_sim)

drug_sim_mat = as.matrix(compound_sim[temp, temp])
##


save(drug_sim_mat, file="kiba_2_drug_sim.rda")

## test
temp1 = dt_triplet[9500,1]
compound_IDs[temp1]
temp2 = dt_triplet[19000,1]
compound_IDs[temp2]
which(df[,1]==compound_IDs[temp1])
which(df[,1]==compound_IDs[temp2])
df[which(df[,1]==compound_IDs[temp1]),2]
df[which(df[,1]==compound_IDs[temp2]),2]
which(compound_ids_sim==df[which(df[,1]==compound_IDs[temp1]),2])
which(compound_ids_sim==df[which(df[,1]==compound_IDs[temp2]),2])
compound_sim[289,1059]
temp3 = which(unique(dt_triplet[,1])==temp1)
temp4 = which(unique(dt_triplet[,1])==temp2)
drug_sim_mat[temp3, temp4]

which(drug_sim_mat[1,]>0.8)
compound_IDs[unique(dt_triplet[,1])[which(drug_sim_mat[1,]>0.8)]]
df[which(df[,1]%in%compound_IDs[unique(dt_triplet[,1])[which(drug_sim_mat[1,]>0.8)]]),2]



