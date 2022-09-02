getCFS <- function(object, ident.1, ident.2, return.object = FALSE){
  CFStbl <- getCFStable(object, ident.1, ident.2)
  # Filter to only contain clusters which have a >40% composition from a single input group
  CFStbl <- CFStbl[, which(apply(CFStbl, MARGIN = 2, FUN = max) > 40)]
  # Assign IDs to new clusters based on majority representation of old cluster
  autoassign <- apply(CFStbl, MARGIN = 2, FUN = function(x) assign.na(x))
  autoassign <- autoassign[!is.na(autoassign)]
  Idents(object) <- plyr::mapvalues(object@meta.data[, ident.2],
                                    from = names(autoassign),
                                    to = rownames(CFStbl)[autoassign])
  # Clusters which were not translated (because < 40% from single group)
  #  get assigned an "NA" value
  low.assign <- which(!levels(object@meta.data[, ident.2]) %in% names(autoassign))
  Idents(object) <- plyr::mapvalues(Idents(object),
                                    from = levels(object@meta.data[, ident.2])[low.assign],
                                    to = rep(NA, length(low.assign))
  )
  # Save translaed IDs as "CFSID"
  object[["CFSID"]] <- Idents(object)
  CFSscore <- sum(as.character(object@meta.data$CFSID) == as.character(object@meta.data[, ident.1]), na.rm = T)/nrow(object[[]])
  if(return.object){
    return(object)
  } else {
    return(round(CFSscore, 4))
  }
}



getCFStable <- function(object, ident.1, ident.2){
  index.df <- data.frame(a = object[[ident.1]], 
                         b = object[[ident.2]])
  index.df <- tidyr::drop_na(index.df)
  colnames(index.df) <- c(paste(ident.1), paste(ident.2))
  index.tbl <- table(index.df)
  # index.tbl <- index.tbl[, which(colSums(index.tbl) > 30)]
  clust <- round(prop.table(index.tbl, margin = 2)*100, 1)
  return(clust)
}


assign.na <- function(x){
  y <- which(x == max(x))
  if(length(y) > 1){
    return(NA)
  } else {
    return(y)
  }
}

