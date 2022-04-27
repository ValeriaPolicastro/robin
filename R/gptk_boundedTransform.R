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



.boundedTransform <-
function (x, transform="atox", bounds) {

  eps <- 2.2204e-16

  thre <- 36	## threshold
  y <- array(0, dim(as.array(x)))

  if ( "atox" == transform ) {
    for ( ind in seq_along(as.array(x)) ) {
      if ( x[ind] > thre )
        y[ind] <- 1-eps
      else if ( x[ind] < -thre )
        y[ind] <- eps
      else
        y[ind] <- 1/(1+exp(-x[ind]))
    }
    y <- (bounds[2] - bounds[1])*y + bounds[1]
  } else if ( "xtoa" == transform ) {
    x <- (x - bounds[1]) / (bounds[2] - bounds[1])
    for ( ind in seq_along(as.array(x)) ) {
      y[ind] <- .complexLog(x[ind]/(1-x[ind]))
    }
  } else if ( "gradfact" == transform ) {
    y <- (x-bounds[1])*(1-(x-bounds[1])/(bounds[2] - bounds[1]))
  }

  return (y)
}
