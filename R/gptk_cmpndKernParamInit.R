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

.cmpndKernParamInit <-
function (kern) {
  
  kern$nParams <- 0
  kern$transforms <- list()

  if ( ! ("comp" %in% names(kern)) )
    kern$comp <- list()

  for ( i in seq(along=kern$comp) ) {

    kern$comp[[i]] <- .kernParamInit(kern$comp[[i]])
    kern$nParams <- kern$nParams + kern$comp[[i]]$nParams
    kern$comp[[i]]$index <- array()

    if ( "numBlocks" %in% names(kern$comp[[i]]) ) {
      if ( i==1 ) {
        kern$numBlocks <- kern$comp[[i]]$numBlocks
      } else {
        if ( (!("numBlocks" %in% names(kern))) | (kern$numBlocks!=kern$comp[[i]]$numBlocks) ) {
          stop("Compound of multi kernels with different numbers of blocks.")
        }
      }
    } else {
      if ( "numBlocks" %in% names(kern) )
        stop("Attempt to combine multi-kernel with non multi-kernel.")
    }
  }

  kern$paramGroups <- diag(1, nrow=kern$nParams, ncol=kern$nParams)

  kern$whiteVariance <- 0
  kern$isStationary <- TRUE

  for ( i in seq(along=kern$comp) ) {
    if ( !kern$comp[[i]]$isStationary )
      kern$isStationary <- FALSE

    if ( kern$comp[[i]]$type == "white" ) {
      kern$whiteVariance <- kern$whiteVariance + kern$comp[[i]]$variance
    } else {
      if ( "whiteVariance" %in% names(kern$comp[[i]]) ) {
        kern$whiteVariance <- kern$whiteVariance + kern$comp[[i]]$whiteVariance
      }
    }
  }

  return (kern)
  
}
