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

#'@importFrom grDevices dev.off dev.set gray.colors pdf rgb
#'@importFrom stats optim optimize qgamma qnorm rnorm
#'@importFrom utils tail
#'@importFrom graphics image

.gpComputeAlpha <-
function(model, m) {

  if (nargs() < 2)
    m = model$m

  model$alpha = matrix(0, model$k, model$d)
  if (model$approx == "ftc") {
# browser()	## m = y-mean(y)
    if (!"isSpherical" %in% names(model) || model$isSpherical)
      model$alpha = model$invK_uu %*% m
    else {
      for (i in 1:model$d) {
	ind = .gpDataIndices(model, i)
	model$alpha[ind, i] = model$invK_uu[[i]] %*% m[ind, i,drop=FALSE]
      }
    }
  }
  else if (model$approx %in% c("dtc","dtcvar")) {
    if (!("isSpherical" %in% names(model)) || model$isSpherical)
      model$alpha = model$Ainv %*% model$K_uf %*% m
    else {
      for (i in 1:model$d) {
	ind = .gpDataIndices(model, i)
	model$alpha[,i] = model$Ainv[[i]] %*% model$K_uf[,ind,drop=FALSE] %*% m[ind,i,drop=FALSE]
      }
    }
  }
  else if (model$approx == "fitc") {
    if (!("isSpherical" %in% names(model)) || model$isSpherical)
      model$alpha = model$Ainv %*% model$K_uf %*% model$Dinv %*% m
    else {
      for (i in 1:model$d) {
	ind = .gpDataIndices(model, i)
	model$alpha[,i] = model$Ainv[[i]] %*%
			    model$K_uf[,ind,drop=FALSE] %*% model$Dinv[[i]] %*% m[ind,i,drop=FALSE]
      }
    }
  } else if (model$approx == "pitc") {
    if (!("isSpherical" %in% names(model)) || model$isSpherical)
      for (i in seq(along=model$blockEnd)) {
	ind = .gpBlockIndices(model, i)
	model$alpha = model$alpha + model$Ainv%*%model$K_uf[,ind,drop=FALSE]%*%model$Dinv[[i]]%*%m[ind, ,drop=FALSE]
      }
    else {
      for (i in seq(along=model$blockEnd)) {
	for (j in 1:model$d) {
	  ind = .gpDataIndices(model, j, i)
	  model$alpha[,j] = model$alpha[,j,drop=FALSE] + model$Ainv[[j]]%*%model$K_uf[,ind,drop=FALSE]%*%
	      model$Dinv[[i]][[j]]%*%m[ind,j,drop=FALSE]
	}
      }
    }
  }

  return(model)
}
