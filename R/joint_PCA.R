library(stats)


# Number of PCA components required to achieve varProp cumulative variance proportion of PCA result full_PCA
cumVar_index <- function(full_PCA, varProp) {
    cumVarProp <- summary(full_PCA)$importance[3,]
    which(cumVarProp > varProp)[1]
}


# Perform joint PCA of training and test data, producing multiple pairs of PCAd datasets according to the
# cumulative variance proportions specified in varProps.
# The options are whether to standardise the data before PCA, whether to plot the components and print the summary
# and whether to include the PCA result object in the output
jointPCA <- function(train, test, varProps = c(0.9), center = T, scale = T, verbose = T, return_full_PCA_data = F) {
    varProps <- sort(varProps)
    names(varProps) <- varProps

    both <- bind_rows(train, test)
    #all.equal(both[1:nrow(train),], train)
    #all.equal(both[(nrow(train) + 1):nrow(both),], test)

    both_PCA <- prcomp(both, center = center, scale. = scale)

    if (verbose) { print(summary(both_PCA)); plot(both_PCA) }

    res <- list(varSets = list())
    if (return_full_PCA_data) { res$full_PCA <- both_PCA }

    train_PCAd <- both_PCA$x[1:nrow(train),]
    test_PCAd <- both_PCA$x[(nrow(train) + 1):nrow(both),]
    
    res$varSets <- map(varProps, function(varProp) {
        ind <- cumVar_index(both_PCA, varProp)
        list(train = train_PCAd[, 1:ind], test = test_PCAd[, 1:ind])
    })
    
    res
}


