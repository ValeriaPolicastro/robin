.gprege <- function(data, inputs, gpregeOptions) {
    
    ## GPREGE Gaussian process ranking and estimation of gene expression time-series.
    ## FORMAT
    ## DESC Fit two GPs with the an RBF (+ noise diagonal) kernel on each
    ## profile. One GP kernel is initialised wih a short lengthscale
    ## hyperparameter, signal variance as the observed variance and a zero noise
    ## variance. It is optimised via scaled conjugate gradients (netlab). The
    ## other GP has fixed hyperparameters with a zero inverse-width, zero signal
    ## variance and noise variance as the observed variance. The log-ratio of
    ## marginal likelihoods of the two hypotheses acts as a score of
    ## differential expression for the profile. Comparison via ROC curves is
    ## performed against BATS (Angelini et.al, 2007).
    ## See Kalaitzis & Lawrence (2011) for a detailed discussion of the
    ## ranking algorithm and dataset used.
    ## ARG data : The matrix of gene expression profiles; one profile per row.
    ## ARG inputs : Inputs to the GP.
    ## ARG gpregeOptions$explore: (LOGICAL) Operate in a user interactive mode.
    ## Used for examining individual gene expression profiles.
    ## ARG gpregeOptions$labels : Contains flags that specify whether the
    ## corresponding profile comes from a differentially expressed gene
    ## (usually from a ground truth).
    ## ARG gpregeOptions$indexRange : Range of indices of profiles on which the
    ## function should operate. Useful for selective exploration of specific
    ## profiles, e.g. only genes marked as differentially expressed in a ground
    ## truth list.
    ## ARG gpregeOptions$interpolatedT : New timepoints to interpolate for each
    ## profile, based on the estimated function values.
    ## ARG gpregeOptions$iters : The number of iterations for scaled-conjugate
    ## gradients (SCG) optimisation.
    ## ARG gpregeOptions$display : Display gradient and LML information on each
    ## SCG iteration.
    ## ARG gpregeOptions$inithypers : The matrix of hyperparameter
    ## configurations as its rows:
    ## [inverse-lengthscale   percent-signal-variance   percent-noise-variance]
    ## The first row corresponds to a (practically constant) function
    ## with a very large lengthscale. Such a function will account for 0 percent
    ## of the observed variance in the expression profile (hence 0 for signal)
    ## and explain it as noise (hence 1 for noise). Subsequent rows
    ## (initialisations for SCG optimisation) correspond to functions of various
    ## lengthscales that explain all the observed variance as signal. A
    ## reasonable lengthscale would be roughly in line with the time-point
    ## sampling intervals.
    ## ARG gpregeOptions$exhaustPlotRes : The search resolution. Used for
    ## interactive mode (explore == 1).
    ## ARG gpregeOptions$exhaustPlotLevels : # Exhaustive plot contour levels.
    ## Used for interactive mode (explore == 1).
    ## ARG gpregeOptions$exhaustPlotMaxWidth : maximum lengthscale to search
    ## for. Used for interactive mode (explore == 1).
    ## RETURN gpregeOutput$signalvar : The vertical lengthscales of the
    ## optimised RBF kernel for each profile.
    ## RETURN gpregeOutput$noisevar : Same, for the noise hyperparameter.
    ## RETURN gpregeOutput$width : Same, for the horizontal lengthscales of the RBF.
    ## RETURN gpregeOutput$LMLs : Log-marginal likelihood of the GP for each profile.
    ## RETURN gpregeOutput$interpolatedData : extended dataset with interpolated values.
    ## RETURN gpregeOutput$rankingScores : the ranking scores based on the
    ## log-ratio of marginal likelihoods.
    ## 
    ## USAGE: gpregeOutput <- gprege(exprs_tp63_RMA, seq(0,240,by=20), gpregeOptions)
    ## 
    ## SEEALSO : gpOptions, gpCreate, gpExpandParam, gpOptimise, gpExtractParam,
    ## gpLogLikelihood, gpPosteriorMeanVar
    ## 
    ## COPYRIGHT: Alfredo A. Kalaitzis, 2010, 2011
    ## 
    ## GPREGE
    
    ## 
    
    if (missing(data))
        stop('Missing data.')
    if (missing(inputs)) {
        stop('Missing inputs.')
    }
    if (is.null(gpregeOptions$indexRange)) {
        gpregeOptions$indexRange <- 1:dim(data)[1]
    }
    n = length(gpregeOptions$indexRange)
    if (is.null(gpregeOptions$explore)) {
        gpregeOptions$explore <- FALSE
    } else {
        if ((gpregeOptions$explore) && is.null(gpregeOptions$exhaustPlotRes)) {
            gpregeOptions$exhaustPlotRes <- 20
        }
        if ((gpregeOptions$explore) && is.null(gpregeOptions$exhaustPlotMaxWidth)) {
            gpregeOptions$exhaustPlotMaxWidth <- 30
        }
    }
    if (is.null(gpregeOptions$iters)) {
        gpregeOptions$iters <- 100
    }
    if (is.null(gpregeOptions$display)) {
        gpregeOptions$display <- FALSE
    }
    gpregeOutput <- list()
    if (!is.null(gpregeOptions$interpolatedT)) {
        interpolate <- TRUE
        newLength = dim(data)[2] + length(gpregeOptions$interpolatedT)
        gpregeOutput$interpolatedData <- matrix(0,newLength, n)
    } else {
        interpolate <- FALSE
    }
    if (is.null(gpregeOptions$inithypers)) {
        gpregeOptions$inithypers <-  matrix( c(1/1000, 0, 1,   1/8, 0.999, 1e-3), ncol=3, byrow=TRUE)
    }
    npsets <- dim(gpregeOptions$inithypers)[1] ## Number of hparams sets.
    if (is.null(gpregeOptions$labels)) {
        gpregeOptions$labels <- rep(0,n)
    }
    
    options <- .gpOptions() ## Set up model.
    options$kern$comp <- list("rbf","white")
    x <- inputs
    xstar <- matrix(seq(min(x)-(2*(max(x)-min(x))/100), max(x)+((max(x)-min(x))/100), length=100), ncol=1)
    
    models <- list() ## Allocate space for vectors.
    loghypers <- matrix(0, 3, npsets)
    LMLs <- matrix(0, n, npsets)
    signalvar <- matrix(0, n, 1); noisevar <- matrix(0, n, 1); width <- matrix(0, n, 1)
    
    ## Remove mean across timepoints. Dataset should not be standardized; Must
    ## look to the signal variance in the context of all signal variances.
    data <- t(as.matrix(scale(t(data), scale=FALSE)))
    
    datamax <- max(data); datamin <- min(data) ## Data min/max for plotting limits.
    
    if (gpregeOptions$explore) {
        #   graphics.off()
        #   dev.new(width=5,height=12);
        #   plot.new() ## Close all devices, open new device, new plot.
        #   dev.new(width=5,height=5);
        #   plot.new();
        #   dev.new(width=5,height=5);
        #   plot.new()
    }
    
    for (ix in 1:n) {
        i <- gpregeOptions$indexRange[ix]
        y <- matrix(data[i,], ncol=1) ## Column vector.
        
        if (sum(is.nan(y)) > (length(y)/2)) {
            cat('Majority of points in profile are NaN.\n')
            next
        }
        
        options$isMissingData <- any(is.nan(y))
        options$isSpherical <- !any(is.nan(y))
        stdy = sd(c(y[!is.nan(y)])) ## Profile variance.
        
        ## Hyperparameters: inverse-lengthscale, signal-variance, noise-variance.
        ## Use 1e-3 instead of 0, or R might coerse them into NaN.
        #   inithypers = t(2 * log(gpregeOptions$inithypers %*% diag(c(1, stdy, stdy))))
        inithypers = t( log( gpregeOptions$inithypers %*% diag(c(1, stdy^2, stdy^2)) ) )
        #   lles = c(1000,8,20)
        #   inithypers = matrix(2*log(c(1e-3, 	1e-3,		 stdy # 1/1.0 (stdy*1e-3) (stdy*0.999)
        # 		      ,	1/8.0,	(stdy*0.999),	(stdy*1e-3)
        # 		      ,	1/20.0,	(stdy*0.999),	(stdy*1e-3)
        # 		      )), nrow=3)
        
        if (gpregeOptions$explore) {
            #     dev.set(2)
            pdf(file=paste('gpPlot',ix,'.pdf',sep=''), paper="special", width=6, height=4*npsets)  
            close.screen(all.screens = TRUE)
            split.screen(c(dim(inithypers)[2],1)) ## Reset any existing sub-figures setup.
            #nf = layout(matrix(c(1:dim(inithypers)[2]), 1:dim(inithypers)[2], 1)); layout.show(nf)
        }
        
        ## Optimise GP log likelihoods.
        for (h in 1:npsets) {
            models[[h]] <- .gpCreate(dim(x)[2], dim(y)[2], x, y, options)
            models[[h]] <- .gpExpandParam(models[[h]], inithypers[, h]) ## This forces kernel computation.
            if (h != 1) {
                models[[h]] <- .gpOptimise(models[[h]], gpregeOptions$display, gpregeOptions$iters)
                loghypers[, h] <- .gpExtractParam(models[[h]], only.values=FALSE)
            }
            LMLs[ix, h] <- .gpLogLikelihood(models[[h]]) ## GP log-marginal likelihood for model h.
            
            if (gpregeOptions$explore) { ## Plot the regression...
                screen(h); erase.screen(h) ## ...in sub-figure #h. Clear screen.
                .gpPlot(models[[h]], xstar, col='blue', xlim=range(xstar), ylim=c(datamin, datamax),
                       title=paste("Init. length-scale = ", as.character(1/exp(inithypers[1,h]/2)), sep=""),
                       xlab='time(mins)', ylab=bquote(paste(log[2]~expression~(centred))))
            }
        }
        
        ## Save maximum log-marginal likelihood and respective hyper-parameters.
        maxLML <- max(LMLs[ix, ]); mi <- which.max(LMLs[ix, ])
        bestLoghypers <- loghypers[, mi]
        width[ix] <- 1/exp(bestLoghypers[1]/2)
        signalvar[ix] <- exp(bestLoghypers[2]/2)
        noisevar[ix] <- exp(bestLoghypers[3]/2)
        
        if (interpolate) {
            meanVar <- .gpPosteriorMeanVar(models[[mi]], gpregeOptions$interpolatedT, varsigma.return=TRUE)
            mu <- meanVar$mu; S <- meanVar$varsigma
            ## Add noise sampled from a Gaussian distribution with zero
            ## mean and variance equal to the predictive variance on new inputs.
            mu <- mu + matrix(rnorm(length(mu)), ncol=1) * noisevar[ix]
            gpregeOutput$interpolatedData[, ix] <- c(y, mu)
            sortedInputs = sort(c(x, gpregeOutput$interpolatedT))
            idx = sortedInputs$ix #sortedInputs = sortedInputs$x
            gpregeOutput$interpolatedData[, ix] = gpregeOutput$interpolatedData[idx, ix] ## Order by augmented inputs.
        }
        
        if (gpregeOptions$explore) {
            if (interpolate) {
                points(gpregeOutput$interpolatedT, mu, pch = 4, cex = .5, lwd=2, col = 'blue')
            }
            cat('\n========================================================\n')
            cat(' Profile ', as.character(i), '\t\t\t\tLabel: ', as.character(as.numeric(gpregeOptions$labels[i])))
            cat('\n========================================================\n')
            cat(sprintf('%-20s %-20s %s\n', 'Length-scale', 'Signal', 'Noise'))
            cat(sprintf('% -20.5f % -20.5f % .5f\n\n', width[ix], signalvar[ix], noisevar[ix]))
            cat(sprintf('%-20s %-20s %s\n', 'Init.le', 'LML', 'Best'))
            best = replace(rep("",npsets), which(LMLs[ix,]==maxLML), " <--")
            for (j in 1:npsets) {
                cat(sprintf('% -20.0f % -20.8f %s\n', 1/exp(inithypers[1,j]/2), LMLs[ix, j], best[j], '\n', sep='\t\t'))
            }
            cat(sprintf('\n%-20s %-20s\n','Total st.dev.','Estim.sig + noise'))
            cat(sprintf('% -20f % -20f\n\n', sd(c(y)), sum(exp(bestLoghypers[2:3]/2))))
            ## We express the profile-ranking metric through a log-ratio of marginal likelihoods.
            cat('Log-ratio (max(LML[2:end]) - LML[1])\n ',max(LMLs[ix,2:npsets])-LMLs[ix,1],'\n')
            
            #     dev.copy(pdf, file = paste('gpPlot',ix,'.pdf', sep=''))
            #     dev.set(3)
            dev.off()
            pdf(file=paste('exhaustivePlot',ix,'.pdf',sep=''), paper="special", width=6, height=6)  
            C = .exhaustivePlot(y, x, xstar, options, gpregeOptions$exhaustPlotMaxWidth, gpregeOptions$exhaustPlotRes, gpregeOptions$exhaustPlotLevels)
            dev.off()
            #     dev.copy(pdf, file = paste('exhaustivePlot',ix,'.pdf', sep=''))
            readline(prompt='ENTER to continue')
        } else{
            cat(' Profile ', as.character(i), '\n')
        }
    }
    
    gpregeOutput$signalvar <- signalvar
    gpregeOutput$noisevar <- noisevar
    gpregeOutput$width <- width
    gpregeOutput$LMLs <- LMLs
    gpregeOutput$rankingScores = apply(as.matrix(LMLs[,2:npsets]),1,max) - LMLs[,1]
    
    return(gpregeOutput)
}