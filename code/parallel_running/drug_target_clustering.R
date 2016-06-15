drug_clusters = hclust(dist(drug_sim_mat))
target_clusters = hclust(dist(sim_mat))

clusterCut_drugs = cutree(drug_clusters,200)
clusterCut_targets = cutree(target_clusters,10)

