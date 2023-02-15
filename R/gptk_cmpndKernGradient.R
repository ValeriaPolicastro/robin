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

.cmpndKernGradient <-
function (kern, x, x2, covGrad) {
  if ( nargs()<4 ) 
    covGrad <- x2
  
  g <- array(0, dim(kern$paramGroups)[1])
  startVal <- 1
  endVal <- 0

  for ( i in seq(along=kern$comp) ) {
    endVal <- endVal + kern$comp[[i]]$nParams
    if ( !is.na(kern$comp[[i]]$index) ) {
      if ( nargs() < 4 ) {
	 g[seq(startVal,len=kern$comp[[i]]$nParams)] <- .kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], covGrad)
        ## Used to be: g[startVal:endVal] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], covGrad)
      } else {
        g[seq(startVal,len=kern$comp[[i]]$nParams)] <- .kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index], covGrad)
        ## Used to be: g[startVal:endVal] <- kernGradient(kern$comp[[i]], x[,kern$comp[[i]]$index], x2[,kern$comp[[i]]$index], covGrad)
      }
    } else {
      if ( nargs() < 4 ) {
        g[seq(startVal,len=kern$comp[[i]]$nParams)] <- .kernGradient(kern$comp[[i]], x, covGrad)
        ## Used to be: g[startVal:endVal] <- kernGradient(kern$comp[[i]], x, covGrad)
      } else {
        g[seq(startVal,len=kern$comp[[i]]$nParams)] <- .kernGradient(kern$comp[[i]], x, x2, covGrad)
        ## Used to be: g[startVal:endVal] <- kernGradient(kern$comp[[i]], x, x2, covGrad)
      }
    }
    startVal <- endVal + 1       
  }

  g <- g %*% kern$paramGroups    

  return (g)
}
