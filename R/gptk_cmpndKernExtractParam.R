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

.cmpndKernExtractParam <-
function (kern, only.values=TRUE,
                                   untransformed.values=FALSE) {

  startVal <- 1
  endVal <- 0
  
  if ( only.values ) {
    params <- c()

    for ( i in seq(along=kern$comp) ) 
      params <- c(params, .kernExtractParam(kern$comp[[i]],
                                           untransformed.values=untransformed.values))

  } else {
    storedTypes <- c()
    params <- c()
    paramNames <- c()
    origNames <- c()
    for ( i in seq(along=kern$comp) ) {
      paramsList <- .kernExtractParam(kern$comp[[i]], only.values=only.values,
                                     untransformed.values=untransformed.values)
      params <- c(params, paramsList)
      kernName <- paste(kern$comp[[i]]$type, length(grep(kern$comp[[i]]$type, storedTypes))+1, sep="")
      paramName <- paste(kernName, names(paramsList), sep="_")
      origNames <- c(origNames, paramName)
      storedTypes <- c(storedTypes, kern$comp[[i]]$type)
    }
  }

  paramNames <- array()
  if ( "paramGroups" %in% names(kern) ) {
    paramGroups <- kern$paramGroups
    for ( i in seq(length.out=dim(paramGroups)[2]) ) {
      ind <- grep(1, paramGroups[,i])
      if ( !only.values ) {
        paramNames[i] <- origNames[ind[1]]
        for ( j in seq(2, length.out=length(ind)-1) )
          paramNames[i] <- paste(paramNames[i], origNames[ind[j]],sep="/")
      }
   
      paramGroups[ind[seq(2,length(ind),length=length(ind)-1)], i] <- 0
    }
  }

  params <- params%*%paramGroups
  if ( !only.values )
    names(params) <- paramNames

  return (params)
}
