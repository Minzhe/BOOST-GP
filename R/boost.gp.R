# ================================== #
#              boost.R               #
# ================================== #

# Load functions
# =========================================================================================
# The core function is boost with six arguments:
#  1) y: an n-by-1 column vector, where n is number of cells/spots and p is the number of genes
#  2) dist: an n-by-n distance matrix, indicating the distance between each spot pair
#  3) nei: a matrix with n rows, where the i-th row lists the neighbors of the i-th spot
#  4) s: a vector with n elements, indicating the size factor for each spot
#  5) kernel: 1 is for squared exponential (SE), 2 is for rational quadratic (RQ), 
#   3 is for Matern, and 4 is for cosine
#  6) model: TRUE for nb model and FALSE for poi model
# The main outputs posterior probability of inclusion (PPI) and original p-value, indicating
#  and the confidence level of being an spatially variable (SV) gene
source("R/functions.R");
Rcpp::sourceCpp("R/boost.cpp");
# =========================================================================================



# function
boost.gp <- function(Y, loc, iter = 1000, burn = 500, size.factor = NA, init.b.sigma = NA, init.h = 1, 
                      update_prop = 0.2, chain = 1, return_lambda = FALSE) {
    res.list <- list()

    for (i in 1:chain) {
        print(paste0("Chain ", i))
        # param
        if (length(size.factor) == 1) {
            if (is.na(size.factor)) {
                size.factor <- rep(1, length(Y))
            }
        }
        if (is.na(init.b.sigma)) {
            init.b.sigma <- round(mean(apply(Y, 2, FUN = function(x) var(log(x[x > 0])))) * 2, 3)
        }
        dist.mat <- distance(loc, 2) # Obtain the distance and neighbor matrix
        
        # run model
        start_time <- proc.time()
        res <- boost(Y = as.matrix(Y), dist = dist.mat$dist, nei = dist.mat$nei, 
                     iter = iter, burn = burn, s = size.factor, 
                     init_b_sigma = init.b.sigma, init_h = init.h, update_prop = update_prop) # See the definition of each argument in the header
        end_time <- proc.time()
        run.time <- as.numeric((end_time - start_time)[3], "secs")
        print(paste0("runtime is ", run.time, "s"))
        
        # significance
        loglr <- res$logBF[(burn+1):iter,]
        loglr[loglr < 0] <- 0
        # p.vals <- apply(as.matrix(loglr), 2, function(x) pchisq(x * 2, df = 1, lower.tail = FALSE))
        # p.vals <- apply(as.matrix(p.vals), 2, function(x) ks.test(x, punif, alternative = "greater")$p.value)
        # p.vals <- combine.pval(p.vals)
        n <- (iter - burn) * update_prop
        bf <- apply(as.matrix(loglr), 2, function(x) mean(sort(x, decreasing = TRUE)[1:n]))
        p.vals <- pchisq(bf * 2, df = 1, lower.tail = FALSE)
        # param
        l <- res$l[(burn+1):iter,]
        # l <- sapply(1:ncol(l), FUN = function(i) mean(l[,i][loglr[,i] > 0]))
        l <- sapply(1:ncol(l), FUN = function(x) sum(l[,x] * loglr[,x]) / sum(loglr[,x]))
        l[is.na(l)] <- 0
        # loglambda
        logLambda <- round(res$logLambda, 3)
        row.names(logLambda) <- paste0("logLambda", 1:nrow(logLambda))
        
        # result
        result <- do.call(rbind, list(l = l, BF = bf, PPI = res$gamma_ppi, pval = p.vals, time = run.time / ncol(Y)))
        colnames(result) <- colnames(Y)
        if (return_lambda) {
            result <- rbind(result, logLambda)
        }
        res.list[[i]] <- data.frame(t(result))
    }
    res.ave <- Reduce("+", res.list) / chain
    pvals <- do.call("rbind", lapply(1:length(res.list), FUN = function(x) res.list[[x]]$pval))
    res.ave$pval <- combine.pval(pvals)
    res.ave$time <- res.ave$time * chain
    return(res.ave)
}
