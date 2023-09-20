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

#' @title TDT Curve plot
#' 
#' 