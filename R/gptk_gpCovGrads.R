.gpCovGrads <-
function(model, M) {

  if (model$approx %in% list('dtc', 'dtcvar')) {
    ## Deterministic training conditional.
    if (model$approx == 'dtcvar')
      dtcvar = TRUE
    else
      dtcvar = FALSE

    if (!'isSpherical' %in% names(model) || model$isSpherical) {
      E = model$K_uf%*%M
      EET = E%*%t(E)
      AinvEET = model$Ainv%*%EET
      AinvEETAinv = AinvEET%*%model$Ainv
      gK_uu = 0.5*(model$d*(model$invK_uu-(1/model$beta)*model$Ainv) - AinvEETAinv)
      if (dtcvar) {
	K_uuInvK_uf = model$invK_uu%*%model$K_uf
	gK_uu = gK_uu - 0.5*model$d*model$beta*(K_uuInvK_uf%*%t(K_uuInvK_uf))
      }
      AinvK_uf = model$Ainv%*%model$K_uf
      gK_uf = -model$d*AinvK_uf - model$beta*(AinvEET%*%AinvK_uf - (model$Ainv%*%E%*%t(M)))
      if (dtcvar)
	gK_uf = gK_uf + model$d*model$beta*K_uuInvK_uf

      gBeta = 0.5*(model$d*((model$N-model$k)/model$beta
		      +sum(model$Ainv * model$K_uu)/(model$beta*model$beta))
		      +sum(AinvEETAinv * model$K_uu)/model$beta
		      +(sum(diag(AinvEET))-sum(M * M)))
      if (dtcvar)
	gBeta = gBeta -0.5*model$d*sum(model$diagD)/model$beta

      fhandle = get(model$betaTransform$func, mode="function")
      gBeta = gBeta*fhandle(model$beta, 'gradfact')
      if (dtcvar)
	g_Lambda = matrix(-0.5*model$beta*model$d, 1, model$N)
      else
	g_Lambda = matrix()
    } else {
      gK_uu = matrix(0, model$k, model$k)
      gK_uf = matrix(0, model$k, model$N)
      gBeta = 0
      for (i in 1:model$d) {
	ind = .gpDataIndices(model, i)
	e = model$K_uf[, ind,drop=FALSE]%*%M[ind, i,drop=FALSE]
	Ainve = model$Ainv[[i]]%*%e
	AinveeT = Ainve%*%t(e)      
	AinveeTAinv = Ainve%*%t(Ainve)
	gK_uu = gK_uu + 0.5*((model$invK_uu - (1/model$beta)*model$Ainv[[i]]) - AinveeTAinv)

	AinvK_uf = model$Ainv[[i]]%*%model$K_uf[, ind,drop=FALSE]
	gK_uf[, ind] = gK_uf[, ind,drop=FALSE] - AinvK_uf
	    - model$beta*(AinveeT%*%AinvK_uf - (Ainve%*%t(M[ind, i, drop=FALSE])))

	gBeta = gBeta + 0.5*(((model$N - model$k)/model$beta
		  +sum(model$Ainv[[i]] * model$K_uu)/(model$beta*model$beta))
		  +sum(AinveeTAinv * model$K_uu)/model$beta
		  +(sum(diag(AinveeT))-sum(M[ind, i,drop=FALSE] * M[ind, i,drop=FALSE])))
      }
      fhandle = get(model$betaTransform$func, mode="function")
      gBeta = gBeta*fhandle(model$beta, 'gradfact')
      g_Lambda = matrix()
    }
  } else if (model$approx == 'fitc') {
    ## Fully independent training conditonal.
    if (!'isSpherical' %in% names(model) || model$isSpherical) {
      E = model$K_uf %*% model$Dinv %*% M
      EET = E%*%t(E)
      AinvE = model$Ainv%*%E
      diagK_fuAinvEMT = t(colSums(model$K_uf * (model$Ainv%*%E%*%t(M))))
      AinvEETAinv = AinvE%*%t(AinvE)
      diagK_ufdAinvplusAinvEETAinvK_fu =
	t(colSums(model$K_uf * ((model$d*model$Ainv + model$beta*AinvEETAinv)%*%model$K_uf)))
      invK_uuK_uf = model$invK_uu%*%model$K_uf
      if (TRUE)
	invK_uuK_ufDinv = invK_uuK_uf%*%model$Dinv
      else
	invK_uuK_ufDinv = solve(t(model$L), model$V)

      diagMMT = rowSums(M * M)
      diagQ = - model$d*model$diagD + model$beta*diagMMT
	      + diagK_ufdAinvplusAinvEETAinvK_fu - 2*model$beta*diagK_fuAinvEMT
      gK_uu = 0.5*(model$d*(model$invK_uu - model$Ainv/model$beta) - AinvEETAinv
	      + model$beta*invK_uuK_ufDinv%*%spam::diag.spam(drop(diagQ))%*%t(invK_uuK_ufDinv))
      gK_uf = -model$beta*invK_uuK_ufDinv%*%spam::diag.spam(drop(diagQ))%*%model$Dinv
	      -model$d*model$Ainv%*%model$K_uf%*%model$Dinv
	      -model$beta*AinvEETAinv%*%model$K_uf%*%model$Dinv
	      +model$beta*model$Ainv%*%E%*%t(M)%*%model$Dinv
      g_Lambda = (0.5*diagQ*model$beta) / (model$diagD * model$diagD)
      gBeta = -sum(g_Lambda)/(model$beta*model$beta)
      fhandle = get(model$betaTransform$func, mode="function")
      gBeta = gBeta*fhandle(model$beta, 'gradfact')
    } else {
      gK_uu = matrix(0, model$k, model$k)
      gK_uf = matrix(0, model$k, model$N)
      g_Lambda = matrix(0, model$N, 1)
      gBeta = 0
      for (i in 1:model$d) {
	ind = .gpDataIndices(model, i)
	K_ufDinvK_uf = model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]]%*%t(model$K_uf[, ind,drop=FALSE])
	e = model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]]%*%M[ind, i,drop=FALSE]
	Ainve = model$Ainv[[i]]%*%e
	AinveeTAinv = Ainve%*%t(Ainve)
	diagK_fuAinveyT = t(colSums(model$K_uf[, ind,drop=FALSE] * (Ainve%*%t(M[ind, i,drop=FALSE]))))
	diagK_ufdAinvplusAinveeTAinvK_fu =
	    t(colSums(model$K_uf[,ind,drop=FALSE] * ((model$Ainv[[i]]
	      + model$beta*AinveeTAinv)%*%model$K_uf[, ind,drop=FALSE])))
	invK_uuK_uf = model$invK_uu%*%model$K_uf[, ind,drop=FALSE]
	invK_uuK_ufDinv = invK_uuK_uf%*%model$Dinv[[i]]
	diagyyT = M[ind, i,drop=FALSE] * M[ind, i,drop=FALSE]
	diagQ = -model$diagD[[i]] + model$beta*diagyyT
		+ diagK_ufdAinvplusAinveeTAinvK_fu - 2*model$beta*diagK_fuAinveyT
	gK_uu = gK_uu + 0.5*(model$invK_uu - model$Ainv[[i]]/model$beta - AinveeTAinv
	  + model$beta*invK_uuK_ufDinv%*%spam::diag.spam(drop(diagQ))%*%t(invK_uuK_ufDinv))
	gK_uf[, ind] = gK_uf[, ind,drop=FALSE]
		-model$beta*invK_uuK_ufDinv%*%spam::diag.spam(drop(diagQ))%*%model$Dinv[[i]]
		-model$Ainv[[i]]%*%model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]]
		-model$beta*AinveeTAinv%*%model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]]
		+model$beta*Ainve%*%t(M[ind, i,drop=FALSE])%*%model$Dinv[[i]]
	g_Lambda[ind] = g_Lambda[ind]
	    + 0.5*model$beta*diagQ / (model$diagD[[i]] * model$diagD[[i]])
      }
      gBeta = gBeta - sum(g_Lambda)/(model$beta*model$beta)
      fhandle = get(model$betaTransform$func, mode="function")
      gBeta = gBeta*fhandle(model$beta, 'gradfact')
    }
  } else if (model$approx == 'pitc') {
    ## Partially independent training conditional.
    if (!'isSpherical' %in% names(model) || model$isSpherical) {
      E = matrix(0, model$k, model$d)
      for (i in 1:length(model$blockEnd)) {
	ind = .gpBlockIndices(model, i)
	E = E + model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]]%*%M[ind, ,drop=FALSE]
      }
      AinvE = model$Ainv%*%E
      AinvEET = AinvE%*%t(E)
      AinvEETAinv = AinvEET%*%model$Ainv
      blockQ = list()
      for (i in 1:length(model$blockEnd)) {
	ind = .gpBlockIndices(model, i)
	K_fuAinvEMT = model$beta*t(model$K_uf[, ind,drop=FALSE])%*%AinvE%*%t(M[ind, ,drop=FALSE])
	blockQ[[i]] = -model$d*model$D[[i]] + model$beta*M[ind, ,drop=FALSE]%*%t(M[ind, ,drop=FALSE])
	    + t(model$K_uf[, ind,drop=FALSE])%*%(model$d*model$Ainv + model$beta*AinvEETAinv)%*%
	    model$K_uf[, ind,drop=FALSE] - K_fuAinvEMT - t(K_fuAinvEMT)
      }
      gK_uu = model$d*model$invK_uu - model$d*model$Ainv/model$beta - AinvEETAinv
      gBeta = 0
      gK_ufBase = -(model$d*model$Ainv + model$beta*AinvEETAinv)%*%model$K_uf
	  + model$beta*AinvE%*%t(M)

      g_Lambda = list()
      gK_uf = matrix(0, dim(gK_ufBase)[1], model$N)
      for (i in 1:length(model$blockEnd)) {
	ind = .gpBlockIndices(model, i)
	invK_uuK_ufDinv = model$invK_uu%*%model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]]
	gK_uu = gK_uu + model$beta*invK_uuK_ufDinv%*%blockQ[[i]]%*%t(invK_uuK_ufDinv)

	gK_uf[, ind] = (gK_ufBase[,ind,drop=FALSE]-model$beta*invK_uuK_ufDinv%*%blockQ[[i]])%*%model$Dinv[[i]]

	g_Lambda[[i]] = 0.5*model$Dinv[[i]]%*%blockQ[[i]]%*%model$Dinv[[i]]*model$beta
	gBeta = gBeta - sum(diag((g_Lambda[[i]]))) / (model$beta*model$beta)
      }
      gK_uu = gK_uu*0.5
      fhandle = get(model$betaTransform$func, mode="function")
      gBeta = gBeta*fhandle(model$beta, 'gradfact')
    } else {
      gK_uu = matrix(0, model$k, model$k)
      gK_uf = matrix(0, model$k, model$N)
      g_Lambda = list()
      for (i in 1:length(model$blockEnd)) {
	if (i == 1)
	  indLen = model$blockEnd[1]
	else
	  indLen = model$blockEnd[i] - model$blockEnd[i-1]

	g_Lambda[[i]] = matrix(0, indLen, indLen)
      }
      gBeta = 0
      for (j in 1:model$d) {
	e = matrix(0, model$k, 1)
	for (i in 1:length(model$blockEnd)) {
	  ind = .gpDataIndices(model, j, i)
	  e = e + model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]][[j]]%*%M[ind, j,drop=FALSE]
	}
	Ainve = model$Ainv[[j]]%*%e
	AinveeT = Ainve%*%t(e)
	AinveeTAinv = AinveeT%*%model$Ainv[[j]]
	blockQ = list()
	for (i in 1:length(model$blockEnd)) {
	  ind = .gpDataIndices(model, j, i)
	  K_fuAinveyT = model$beta*t(model$K_uf[, ind,drop=FALSE])%*%Ainve%*%t(M[ind, j,drop=FALSE])
	  blockQ[[i]] = -model$D[[i]][[j]] + model$beta*M[ind, j,drop=FALSE]%*%t(M[ind, j,drop=FALSE])
	      + t(model$K_uf[, ind,drop=FALSE])%*%(model$Ainv[[j]] + model$beta*AinveeTAinv)%*%
		model$K_uf[, ind,drop=FALSE] - K_fuAinveyT - t(K_fuAinveyT)
	}
	gK_uu = gK_uu + model$invK_uu - model$Ainv[[j]]/model$beta - AinveeTAinv

	for (i in 1:length(model$blockEnd)) {
	  ind = .gpDataIndices(model, j, i)
	  gK_ufBase = -(model$Ainv[[i]] + model$beta*AinveeTAinv)%*%model$K_uf[, ind,drop=FALSE]
	      + model$beta*Ainve%*%t(M[ind, j,drop=FALSE])
	  invK_uuK_ufDinv = model$invK_uu%*%model$K_uf[, ind,drop=FALSE]%*%model$Dinv[[i]][[j]]
	  gK_uu = gK_uu + model$beta*invK_uuK_ufDinv%*%blockQ[[i]]%*%t(invK_uuK_ufDinv)

	  gK_uf[, ind] = gK_uf[, ind,drop=FALSE] + (gK_ufBase
			  -model$beta*invK_uuK_ufDinv%*%blockQ[[i]])%*%model$Dinv[[i]][[j]]

	  if (i == 1)
	    localInd = ind
	  else
	    localInd = ind - (model$blockEnd[i-1])

	  g_Lambda[[i]][localInd, localInd] = g_Lambda[[i]][localInd, localInd]
	      + 0.5*model$Dinv[[i]][[j]]%*%blockQ[[i]]%*%model$Dinv[[i]][[j]]*model$beta
	}
      }

      for (i in 1:length(model$blockEnd))
	gBeta = gBeta - sum(diag((g_Lambda[[i]])))/(model$beta*model$beta)

      gK_uu = gK_uu*0.5
      fhandle = get(model$betaTransform$func, mode="function")
      gBeta = gBeta*fhandle(model$beta, 'gradfact')
    }
  } else
    stop("Unknown approximation type")

  return (list(gK_uu=gK_uu, gK_uf=gK_uf, g_Lambda=g_Lambda, gBeta=gBeta))
}
