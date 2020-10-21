library(MASS)
library(ggplot2)

quality_controller <- function(count, loc, cutoff_sample = 10, cutoff_feature = 0.1) {
    n <- dim(count)[1]
    p <- dim(count)[2]
    index <- which(rowSums(count) >= cutoff_sample)
    count <- count[index,]
    loc <- loc[index,]
    index <- which(colSums(count != 0) >= dim(count)[1]*cutoff_feature)
    count <- count[, index]
    return (list("count" = count, "loc" = loc))
}

distance <- function(loc, cutoff) {
    dist <- as.matrix(dist(loc, method = "euclidean", diag = TRUE, upper = TRUE))
    rownames(dist) <- NULL
    colnames(dist) <- NULL
    
    n <- dim(dist)[1]
    neighbors <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
        temp <- c(which(dist[i,] <= cutoff) - 1, -1)
        neighbors[i, 1:length(temp)] <- temp
    }
    neighbors <- neighbors[, which(colSums(neighbors) != 0)]
    return (list("dist" = dist, "nei" = neighbors))
}

gp_generator <- function(dist, mean = rep(0, dim(dist)[1]), sigma = 1, kernel = c("C", "SE", "RQ", "MA", "PE", "COS"), c = 0, l = 1, alpha = 1, p = 1, psi = 1, seed = 1, y_known = rep(NA, dim(dist)[1])) {
    n <- dim(dist)[1];
    if (kernel == "C") {
        if (c >= 1 || c < 0) {
            stop("Invalid value of c!")
        } else {
            K <- matrix(c, nrow = n, ncol = n);
            diag(K) <- 1;
        }
    } else if (kernel == "SE") {
        K <- exp(-(dist^2)/2/l^2);
    } else if (kernel == "RQ") {
        K <- (1 + dist^2/2/alpha/l^2)^(-alpha);
    } else if (kernel == "MA") {
        temp <- 0;
        nu <- p + 1/2;
        for (i in 0:p) {
            temp <- temp + factorial(p + i)/factorial(i)/factorial(p - i)*(sqrt(8*nu)*dist/l)^(p - i);
        }
        K <- exp(-sqrt(2*nu)*dist/l)*gamma(p + 1)/gamma(2*p + 1)*temp;
    } else if (kernel == "PE") {
        K <- exp(-2*sin(pi*dist/psi)^2/l^2);
        temp <- eigen(K);
        P <- temp$vectors;
        D <- diag(temp$values);
        diag(D)[which(diag(D) < 0)] <- 0;
        K <- P %*% D %*% solve(P);
    } else if (kernel == "COS") {
        K <- cos(2*pi*dist/psi);
        temp <- eigen(K);
        P <- temp$vectors;
        D <- diag(temp$values);
        diag(D)[which(diag(D) < 0)] <- 0;
        K <- P %*% D %*% solve(P);
    } else {
        stop("Invalid kernel name!");
    }
    
    set.seed(seed);
    if (sum(is.na(y_known)) == n) {
        y <- mvrnorm(1, mean, Sigma = sigma^2*K);
    } else {
        y <- rep(NA, n);
        index_known <- which(!is.na(y_known));
        y[index_known] <- y_known[index_known];
        mean_unknown <- mean[-index_known] + K[-index_known, index_known] %*% solve(K[index_known, index_known]) %*% (y_known[index_known] - mean[index_known]);
        K_unknown <- K[-index_known, -index_known] - K[-index_known, index_known] %*% solve(K[index_known, index_known]) %*% K[index_known, -index_known]
        y[-index_known] <- mvrnorm(1, mean_unknown, Sigma = sigma^2*K_unknown);
    }
    return (list("y" = y, "K" = K, "mean" = mean, "sigma" = sigma));
}

anscombe_transformer <- function(count) {
    var <- apply(count, 2, var)
    mean <- apply(count, 2, mean)
    phi <- (sum(mean^2*(var - mean)))/sum(mean^4)
    count <- log(count + 1/phi)
    return (count)
}

normalize <- function(count, mode = c("linear", "sigmoid"), alpha = 5) {
    if (mode == "linear") {
        max <- max(count)
        min <- min(count)
        count <- (count - min)/(max - min)
    } else if (mode == "sigmoid") {
        median <- median(count)
        count <- count - median
        count <- 1/(1 + exp(-count*alpha))
    }
    return (count)
}

loglklh <- function(y, dist, l, X = matrix(1, nrow = length(y), ncol = 1), h = 1000, alpha_sigma = 3, beta_sigma = 1) {
    p <- dim(X)[1]
    R <- dim(X)[2]
    if (l == 0) {
        H <- diag(p)
    } else {
        H <- exp(-(dist^2)/2/l^2)
    }
    res <- -Inf
    
    tryCatch({
        Hi <- chol2inv(chol(H))
        # Hi <- solve(H); 
        if (R == 1) {
            G <- 1/h + sum(Hi)
            Gi <- 1/G
            F <- Hi - (Gi)*(Hi %*% matrix(1, nrow = p, ncol = p) %*% Hi)
        } else {
            G <- diag(R)/h + (t(X) %*% Hi %*% X)
            Gi <- chol2inv(chol(G))
            # Gi <- solve(G)
            F <- Hi - (Hi %*% X %*% Gi %*% t(X) %*% Hi)
        }
        res <- -determinant(H, logarithm = TRUE)$modulus/2 - (alpha_sigma + p/2)*log(beta_sigma + (t(y) %*% F %*% y)/2 )
        if (R == 1) {
            res <- res - log(G)/2
        } else {
            res <- res - determinant(G, logarithm = TRUE)$modulus/2
        }
    }, error = function(e){})
    
    return(as.numeric(res))
}

combine.pval <- function(pvals, weights = NULL){
    if(!is.matrix(pvals)) {
        pvals <- as.matrix(pvals)
    }
    # avoid extremely values
    pvals[which(pvals == 0)] <- 5.55e-17
    pvals[which((1 - pvals) < 1e-3)] <- 0.99
    
    n.pval <- nrow(pvals)
    n.gene <- ncol(pvals)
    if(is.null(weights)) {
        weights <- matrix(rep(1.0 / n.pval, n.pval * n.gene), ncol = n.gene )
    }
    if((nrow(weights) != n.pval) || (ncol(weights) != n.gene)){
        stop("The dimensions of weights does not match that of combined pvalues")
    }
    
    Cstat <- tan((0.5 - pvals) * pi)
    wCstat <- weights * Cstat
    Cbar <- apply(wCstat, 2, sum)
    #combined.pval <- 1.0/2.0 - atan(Cbar)/pi
    combined.pval <- 1.0 - pcauchy(Cbar)	
    combined.pval[which(combined.pval <= 0)] <- 5.55e-17
    return(combined.pval)
}

plot.expr <- function(count, loc, main = "") {
    data <- data.frame(expr = as.vector(count), x = loc[, 1], y = loc[, 2]);
    ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = expr), size = 4) + 
        coord_fixed(ratio = 1) + scale_color_distiller(palette = "Spectral") + 
        theme_classic() + labs(color = "", title = main) + 
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              plot.title=element_text(hjust = 0.5),
              legend.position="right")
}
