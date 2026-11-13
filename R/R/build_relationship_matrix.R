# 09/24

library(dplyr)
library(AGHmatrix)


# for later
#' @description builds additive genomic relationship matrix from marker data. Pipeline calls a prebuilt version from `data/processed/markers/G_mat.rds` 
#' @param snp_matrix Genotype x marker snp dosage matrix
#' @param ploidy 4 for autotetraploid
#' @param method - default is `VanRaden`
#' @param missing_value - code used for missing marker calls
#' @param maf minor allele freq filtering threshold
#' @output the additive genomic relationship matrix, dim is genotype x genotype 
build_genomic_g_matrix = function(snp_matrix, ploidy = 4, method = 'VanRaden', missing_value = NA, maf = 0.05) {
  Gmatrix(SNPmatrix = snp_matrix, method = method, ploidy = ploidy, missingValue = missing_value, maf = maf)
}


#' @description builds the pedigree-based additive relationship matrix
#' @param ped_dta a processed pedigree table with three columns indicating, name, male, and female ids
#' @param ploidy defaults to 4
#' @param w inbreeding/weighting option. defaults to 0
#' @output the pedigree based additive relationship matrix, dim is genotype x genotype 
build_pedigree_a_matrix = function(ped_dta, ploidy = 4, w = 0){

  ped_dta = ped_dta %>% rename(Ind = NAME, Par1 = FEMALE, Par2 = MALE)
  Amatrix(data = as.data.frame(ped_dta), ploidy = ploidy)
}


# 10/24 -- genomic dominance relationship matrix 
#' @description builds the dominance genomic relationship matrix from marker data with the Endelman dominance method. Pipeline actually reads a pre-built version from `data/processed/markers/D_mat.rds`
#' @param snp_matrix Genotype x marker snp dosage matrix
#' @param ploidy 4 for autotetraploid
#' @param method - default is `Endelman`
#' @param missing_value - code used for missing marker calls
#' @param maf minor allele freq filtering threshold
#' @output the dominance genomic relationship matrix, dim is genotype x genotype 
build_genomic_d_matrix = function(snp_matrix, ploidy = 4, method = 'Endelman', missing_value = NA, maf = 0.05) {
  Gmatrix(SNPmatrix = snp_matrix, method = method, ploidy = ploidy, missingValue = missing_value, maf = maf)
}


# standard diagonal bending for g_matrix that could not invert reliably
#' @description diagonal bending of a relationship matrix, shrinks its own diagonal by weight `w`
#' @param G relationship matrix, dim genotype x genotype
#' @param w bending weight [0, 1] - larger values pull the matrrix more strongly to its diagonal
#' @output the bent relationship matrix, same dim as G
bend_G = function(G, w = 1e-3) {
  (1 - w) * G + w * diag(diag(G))
}

# mat is A or G mat
#' @description subsets a relationship matrix before inversing
#' @param mat: output of `build_genomic_g_matrix()` or `build_pedigree_a_matrix()`
#' @param keep_genotypes character vector of genotype names to retain for the model.
#' @output the inverse of `mat` with only the genotypes in `keep_genotypes`
build_A_inverse = function(mat, keep_genotypes){
  stopifnot(!is.null(rownames(mat)), all(keep_genotypes %in% rownames(mat)))
  keep_mat = mat[keep_genotypes, keep_genotypes]

  solve(keep_mat)
}


#' @description pulls out the diagonal of rel mat for a given set of genos
#' @param mat: output of `build_genomic_g_matrix()` or `build_pedigree_a_matrix()`
#' @param keep_genotypes character vector of genotype names to retain for the model.
#' @output named numeric vector of the mat diagonal values for `keep_genotypes`
get_a_diag = function(mat, keep_genotypes){
  stopifnot(!is.null(rownames(mat)), all(keep_genotypes %in% rownames(mat)))
  
  setNames(diag(mat), rownames(mat))[keep_genotypes]
}