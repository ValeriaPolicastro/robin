.exhaustivePlot <-
function(y, x, xstar, options, maxwidth, res, nlevels) {

## EXHAUSTIVEPLOT Plot of the LML function by exhaustive search.
## FORMAT 
## DESC Exhaustively searches the hyperparameter space by a grid, whose
## resolution is given, and plots the LML function for every point in the
## space.
## ARG y : the target (output) data.
## ARG x : the input data matrix.
## ARG xstar : the points to predict function values.
## ARG options : options structure as defined by gpOptions.m.
## ARG maxWidth : maximum lengthscale to search for.
## ARG res : The search resolution. Number of points to plot for in the
## search range.
## ARG nlevels : Number of contour levels.
## RETURN C : Matrix of function values from the search.
## 
## USAGE : exhaustivePlot(y, x, xstar, options, 30, 100);
## 
## SEEALSO : gpCreate, gpExpandParam, gpLogLikelihood, gpPosteriorMeanVar
## 
## COPYRIGHT : Alfredo Kalaitzis, 2010, 2011
## 
## GPREGE
  
y = y[!is.nan(y)]
y = matrix(y-mean(y))

y_sd = sd(c(y))
if (y_sd == 0) {
  C <- NULL
  warning('Data variance is zero. No figure produced.')
} else {
  model = .gpCreate(dim(x)[2], dim(y)[2], x, y, options)

  ## search GP log likelihood
  signal = seq(y_sd*1e-3, y_sd*.999, length=res)
  width = seq(1, maxwidth, length=res)
  results = matrix(0, length(signal)*length(width), 5)

  index = 0
  for (w in width) {
    for (sig in signal) {
      noise = y_sd - sig
      snr = sig/noise
      model = .gpExpandParam(model, matrix(2*log(c(1/w, sig, noise)), ncol=1) )
      LML = .gpLogLikelihood(model)
      index = index + 1
      results[index, ] = c(w, sig, noise, snr, LML)
    }
  }

  C = matrix(results[, 5], length(signal), length(width))
  # lbound = max(C) - ((max(C)-min(C))/10)
  lbound = -20;
  C[C < lbound] = lbound ## Put a lower bound on the LML plot.
  maxLML = max(results[,5]); mi = which.max(results[,5])
  v = seq(min(C)+10, maxLML, length=nlevels)
  w = results[mi,1]; sig = results[mi,2]; noise = results[mi,3]; snr = results[mi,4]
  # screen(1); erase.screen(1) ## Subfigure 1; Clear screen.
  ## Plot contour of log-marginal likelihood function wrt s/n ratio.
  ## NOTE: filled.contour maps the transpose of the matrix to the image, hence the t(C).
  filled.contour(x=width, y=log10(signal/(y_sd-signal)), z=t(C),
  # filled.contour(x=width, y=(signal/(sd(y)-signal)), z=C,
    color.palette=gray.colors, #terrain.colors,
      xlab='Lengthscale', ylab=bquote(paste(log[10]~SNR)),
  #   levels = v,
    nlevels = nlevels,
  #   plot.title=title(main='Log-marginal likelihood function\n lengthscale vs. log10(SNR)'),
    main='Log-marginal likelihood function',
    plot.axes={axis(1); axis(2);
      points(w, log10(snr), pch=20)}
  )


  ## plot GP regression with maxLML hyperparameters
  loghyper = matrix(2*log(c(1/w, sig, noise)), ncol=1)
  model = .gpExpandParam(model, loghyper)
  dev.set(4) # screen(2); erase.screen(2)
  .gpPlot(model, xstar, col='blue', # xlim=range(xstar), ylim=c(min(y), max(y)),
    title=paste("max LML length-scale = ", as.character(round(w,1)), sep=''),
    xlab='time(mins)', ylab=bquote(paste(log[2]~expression~(centred))))

  ## Hyperparameter info.
  cat('============= EXHAUSTIVE LML SEARCH =====================\n')
  cat( sprintf('%-20s %-20s %s\n', 'Length-scale', 'Signal', 'Noise') )
  cat( sprintf('%-20.5f %-20.5f %.5f\n\n', w, sig, noise) )
  cat( sprintf('%-20s %-20s\n', 'Max LML', 'Estim. sig + noise') )
  cat( sprintf('%-20.8f %-20f\n\n', gpLogLikelihood(model), sum(exp(loghyper[2:3]/2))) )
}
return(C)
}
