#Tree checks function
tree_checks <- function(data, tree, dataCol, type = c("checks", "prune")){
  type = match.arg(type)
  # How many unique species exist in data and tree
  Numbers <- matrix(nrow = 2, ncol = 1)
  Numbers[1,1] <- length(unique(data[,dataCol]))
  Numbers[2,1] <- length(tree$tip.label)
  rownames(Numbers)<- c("Species in data:", "Species in tree:")
  # Missing species or species not spelt correct
  species_list1= setdiff(sort(tree$tip.label), sort(unique(data[,dataCol])))
  species_list2= setdiff(sort(unique(data[,dataCol])), sort(tree$tip.label) )
  if(type == "checks"){
    return(list(SpeciesNumbers = data.frame(Numbers),
                Species_InTree_But_NotData=species_list1,
                Species_InData_But_NotTree=species_list2))
  }
  if(type == "prune"){
    if(length(species_list2) >=1) stop("Sorry, you can only prune a tree when you have no taxa existing in the data that are not in the tree")
    return(ape::drop.tip(tree, species_list1))
  }
}


#' @title TDT Curve Calculation
#' @description This function calculates the slope and CTmax from a TDT curve fit using a linear regression
#' @param data A dataframe containing the data  
#' @return A data frame containing the slope and CTmax
#' @export 
tdt_calc <- function(data){

  # Fit a linear regression to the data
      fit <- lm(log10(t_coma) ~ assay_temp, data = data)
  
  # Extract the slope and CTmax
   slope <- -1/fit$coefficients[2]
   CTmax <- -coef(fit)[1] / coef(fit)[2]
  
  # Return the slope and CTmax
    return(data.frame(z = slope, CTmax = CTmax))

}

#' @title extract_blups
#' @description This function extracts the BLUPs from a mixed model
#' @param sol A matrix containing the BLUPs
#' @param trait A string. The trait to be analysed. 
#' @return A list containing the phylogenetic BLUPs and the species BLUPs
extract_blups <- function(sol, trait){
  animal <- paste0("trait", trait, ".animal")
  species <- paste0("trait", trait, ".species")

  phy <- sol[, grep(animal, colnames(sol))]
  phy <- phy[,-grep("Node", colnames(phy))]

  spp <- sol[, grep(species, colnames(sol))]

  return(list(phy = phy, spp = spp))

}
  

#' @title predict_accuracy
#' @description This function calculates the accuracy of the model predictions compared to the data where we know the actual values
#' @param data A dataframe containing the data
#' @param blups A list containing the phylogenetic BLUPs and the species BLUPs
#' @param trait The trait to be analysed. Note that this is not a string, but the name of the trait
#' @return A data frame containing the actual and predicted values for a given trait
predict_accuracy <- function(data, blups, trait){
  # Get actual data. Note that we have species replicated multiple times. Stick to among species level so need to average the replicate populations of a given species. BLUPs are only estimated for a single species. 

  known <- data %>% filter(!is.na({{trait}})) %>% dplyr::select(species, {{trait}}) %>% group_by(species) %>% summarise(trait = mean({{trait}}))
  
  # Get estimated data
    est <- apply(blups, 2, function(x) median(x))
    est <- est[which(colnames(blups) %in% known$species)]
    est <- data.frame(est)
    est$species <- row.names(est)

  # Join the two
  check <- left_join(known, est, by = "species")

  return(check)
}

#' @title phylo_corrs
#' @description This function calculates the phylogenetic correlation between two traits
#' @param vcv A variance covariance matrix
#' @param trait1 The first trait to be analysed
#' @param trait2 The second trait to be analysed
#' @return A vector containing the posterior distribution of the phylogenetic correlation between the two traits
phylo_corrs <- function(vcv, trait1, trait2){
    
    # Grep strings
    animalt1 <- paste0("trait", trait1)
    animalt2 <- paste0("trait", trait2)
    
    # Extract the variance components
    vart1 <- vcv[,grep(paste0(animalt1, ":", animalt1, ".animal"), colnames(vcv))]
    vart2 <- vcv[,grep(paste0(animalt2, ":", animalt2, ".animal"), colnames(vcv))]

    # extract covariance between traits
    cov_name <- paste0(animalt1, ":", animalt2, ".animal")

    # Calculate correlatoon
    r_phylo <- vcv[,cov_name] / sqrt(vart1 * vart2)

    return(r_phylo)
}