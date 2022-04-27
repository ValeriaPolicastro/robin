# Based on http://opensource.org/licenses/BSD-2-Clause
#
# YEAR: 2013
# COPYRIGHT HOLDER: Alfredo Kalaitzis
#
#     Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#     
#     Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the
# distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

.CGoptim <-
function (x, fn, grad, options, ...) {
  ## option[1] : number of iterations
  ## option[2] : interval for the line search
  ## option[3] : tolerence for x to terminate the loop
  ## option[4] : tolerence for fn to terminate the loop
  ## option$display : option of showing the details of optimisaton

  ## y = fn (x)
  func <- function(x, ...) fn(x, ...)
      
  ## gradient function = gr (x)
  gradfunc <- function(x, ...) grad(x, ...)
	
  fn_new <- func(x, ...)
  #if ( display ) 
  #  cat ("fn0 :",fn_new, "\n")

  grad_new <- gradfunc(x, ...)
  #if ( options$display )  
  #  cat ("grad0 :",grad_new, "\n\n")
	
  direction <- -grad_new
  lnSchFail <- FALSE
  for ( ind in 1:options[[1]] ) {
	
    x_old <- x
    fn_old <- fn_new
    grad_old <- grad_new
    
    grad2 <- crossprod(grad_old)	
    
    if ( grad2 == 0 ) {
      objective <- fn_new
      xmin <- x
      ans <- list(xmin=xmin, objective=objective, lnSchFail=lnSchFail)
      return (ans)
    }

    dnorm <- sqrt(sum(direction*direction))
    line_dir <- direction / dnorm
    ## cat ("\n line_dir :", line_dir, "\n\n")
    lnSch <- try( optimize(.fn_line, options[[2]], para0=x_old, direction=line_dir, fun=fn, ...) )

    if ( is.list(lnSch) ) {
      x <- x_old + lnSch$minimum * line_dir		
      fn_new <- lnSch$objective
      fnmin <- min(fn_old, fn_new)
      if ( fnmin==fn_new ) {
        xnmin <- x
      } else {
        xnmin <- x_old
      }
    } else {
      warning("Line search failed! \n")
      x <- xnmin
      fn_new <- fnmin
      lnSchFail <- TRUE

      xmin <- xnmin
      objective <- fnmin
      ans <- list(xmin=xmin, objective=objective, lnSchFail=lnSchFail)
      return (ans)
    }
    
    if ( max(abs(x-x_old))<options[[3]] & max(abs(fn_new-fn_old))<options[[4]] ) {
      xmin <- x
      objective <- fn_new
      ans <- list(xmin=xmin, objective=objective, lnSchFail=lnSchFail)
      return (ans)
    }
    
    grad_new <- gradfunc(x, ...)
    
    eta <- ( t(grad_new-grad_old) %*% grad_new ) / grad2
    direction <- direction * eta - grad_new

    if ( options$display )
      cat(ind, "-th objective = :", fn_new, "\t max xi: ", max(abs(x-x_old)), "\n")
  }

  warning("Maximum iteration reached! \n")
  xmin <- x
  objective <- fn_new
  
  ans <- list(xmin=xmin, objective=objective, lnSchFail=lnSchFail)
  ##if ( options$display ) 
  ##  cat("\n Optimisation ends! \n")
  return (ans)
}
