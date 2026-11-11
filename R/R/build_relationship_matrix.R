# 09/24

library(dplyr)
library(AGHmatrix)


# for later
build_genomic_g_matrix = function(snp_matrix, ploidy = 4, method = 'VanRaden', missing_value = NA, maf = 0.01) {
  Gmatrix(SNPmatrix = snp_matrix, method = method, ploidy = ploidy, missingValue = missing_value, maf = maf)
}


build_pedigree_a_matrix = function(ped_dta, ploidy = 4, w = 0){

  ped_dta = ped_dta %>% rename(Ind = NAME, Par1 = FEMALE, Par2 = MALE)
  Amatrix(data = as.data.frame(ped_dta), ploidy = ploidy)
}


# 10/24 -- genomic dominance relationship matrix 
build_genomic_d_matrix = function(snp_matrix, ploidy = 4, method = 'Endelman', missing_value = NA, maf = 0.01) {
  Gmatrix(SNPmatrix = snp_matrix, method = method, ploidy = ploidy, missingValue = missing_value, maf = maf)
}


# standard diagonal bending for g_matrix that could not invert reliably
bend_G = function(G, w = 1e-3) {
  (1 - w) * G + w * diag(diag(G))
}

# mat is A or G mat
#'@param mat: output of build_genomic_g_matrix or build_pedigree_a_matrix
build_A_inverse = function(mat, keep_genotypes){
  stopifnot(!is.null(rownames(mat)), all(keep_genotypes %in% rownames(mat)))
  keep_mat = mat[keep_genotypes, keep_genotypes]

  solve(keep_mat)
}

get_a_diag = function(mat, keep_genotypes){
  stopifnot(!is.null(rownames(mat)), all(keep_genotypes %in% rownames(mat)))
  
  setNames(diag(mat), rownames(mat))[keep_genotypes]
}