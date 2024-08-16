
summarize.subgroups <- function(x, ...) UseMethod("summarize.subgroups")
summarize.subgroups.default <- function(x, subgroup, ...)
{
    vnames <- colnames(x)

    n.obs  <- NROW(x)
    n.vars <- NCOL(x)

    if (is.null(vnames))
    {
        vnames <- paste0("V", 1:n.vars)
    }


    # find which variables are binary
    var.levels <- numeric(n.vars)
    for (v in 1:n.vars)
    {
        var.levels[v] <- length(unique(x[,v]))
    }

    contin.vars <- vnames[var.levels > 2]
    binary.vars <- vnames[var.levels == 2]


    unique.trts <- sort(unique(subgroup))
    n.trts      <- length(unique.trts)

    if (n.trts < 2) stop("There is only one unique subgroup. No subgroups to compare with.")

    compare.mat <- array(0, dim = c(n.vars, 2 * (n.trts + choose(n.trts, 2))))
    colnames(compare.mat) <- 1:ncol(compare.mat)

    for (t in 1:n.trts)
    {
        ## means within each subgroup
        compare.mat[,t] <- colMeans(x[subgroup == unique.trts[t], ])
    }

    for (v in 1:n.vars)
    {

        ct <- 0
        for (t in 1:n.trts)
        {

            if (var.levels[v] > 2)
            {
                ## standard errors within each subgroup
                compare.mat[v,n.trts + 2 * choose(n.trts, 2) + t] <-
                    sd(x[subgroup == unique.trts[t], v]) / sqrt(sum(subgroup == unique.trts[t]))

                if (t < n.trts)
                {
                    for (k in (t + 1):n.trts)
                    {
                        ct <- ct + 1
                        ## run t.test for contin vars
                        tt <- t.test(x[subgroup == unique.trts[t], v], x[subgroup == unique.trts[k], v])
                        compare.mat[v, n.trts + choose(n.trts, 2) + ct] <- tt$p.value
                    }
                }
            } else
            {
                if (t < n.trts)
                {
                    for (k in (t + 1):n.trts)
                    {
                        ct <- ct + 1

                        sub.idx <- subgroup == unique.trts[t] | subgroup == unique.trts[k]
                        ## run chi squared test for binary vars
                        if (length(unique(x[sub.idx, v])) > 1 & sum(sub.idx) > 2)
                        { 
                          cst <- chisq.test(subgroup[sub.idx], x[sub.idx, v])
                          compare.mat[v, n.trts + choose(n.trts, 2) + ct] <- cst$p.value
                            #messg <- tryCatch(cst <- chisq.test(subgroup[sub.idx], x[sub.idx, v]),
                                              #warning = function(w) return(w))
                          
# 
#                             if(is(messg[[2]], "warning"))
#                             {
#                                 #cst <- chisq.test(subgroup[sub.idx], x[sub.idx, v],
#                                 #                  simulate.p.value = TRUE)
#                                 cst <- fisher.test(subgroup[sub.idx], x[sub.idx, v])
#                             }

                            
                        } else
                        {
                            compare.mat[v, n.trts + choose(n.trts, 2) + ct] <- NA
                        }
                    }
                }
            }

        }

    }

    ct <- 0
    for (t in 1:n.trts)
    {
        colnames(compare.mat)[t] <- paste0("Avg (recom ", unique.trts[t], ")")
        colnames(compare.mat)[n.trts + 2 * choose(n.trts, 2) + t] <- paste0("SE (recom ", unique.trts[t], ")")
        if (t < n.trts)
        {
            for (k in (t + 1):n.trts)
            {
                ct <- ct + 1
                compare.mat[,n.trts + ct] <- compare.mat[,t] - compare.mat[,k]
                colnames(compare.mat)[n.trts + ct] <- paste0(unique.trts[t], " - ", unique.trts[k])
                colnames(compare.mat)[n.trts + choose(n.trts,2) + ct] <-
                    paste0("pval ", unique.trts[t], " - ", unique.trts[k])

                compare.mat[,n.trts + choose(n.trts,2) + ct] <- stats::p.adjust(compare.mat[,n.trts + choose(n.trts,2) + ct],
                                                                                "BH")
            }
        }
    }

    rownames(compare.mat) <- vnames

    compare.mat <- as.data.frame(compare.mat)
    #colnames(compare.mat) <- c("avg (recom trt)", "avg (recom ctrl)", "diff",
    #                           "p.value", "SE (recom trt)", "SE (recom ctrl)")
    class(compare.mat)    <- c("subgroup_summary", "data.frame")
    compare.mat
}


summarize.subgroups.subgroup_fitted <- function(x, ...)
{

    if (is.null(x$call)) stop("retcall argument must be set to TRUE for fitted model
                                    to use summarize.subgroups()")


    # save data objects because they
    # will be written over by resampled versions later
    xx       <- x$call$x
    subgroup <- x$recommended.trts

    vnames   <- x$var.names

    colnames(xx) <- vnames
    summarize.subgroups.default(x = xx, subgroup = subgroup)
}


print.subgroup_summary <- function(x, p.value = 0.001, digits = max(getOption('digits')-3, 3), ...)
{
    pidx <- grep("pval", colnames(x))
    lessthan <- x[,pidx,drop = FALSE] <= p.value
    lessthan[is.na(lessthan)] <- FALSE
    if (!is.null(dim(lessthan)))
    {
        compare.mat <- x[rowSums(lessthan) > 0,]
    } else
    {
        compare.mat <- x[lessthan > 0,]
    }
    print.data.frame(compare.mat[,-pidx], digits = digits, quote = FALSE, right = TRUE, na.print = "NA", ...)
}

