#possible that lm model returns NAs in coefficients
#makes diagnostics incalculable
fastVARX = function(y, x, p, b, getdiag=T) {
  if(p < 1) stop("p must be a positive integer")
  if(missing(x)) {
    return (fastVAR(y, p, getdiag))
  }
  varxz = VARXZ(y, x, p, b)
  model = lm(varxz$y.p ~ varxz$Z)
  if(sum(is.na(model$coefficients)) > 0) {
    stop("Multivariate lm has invalid coefficients.  
      Check the rank of the design matrix")
  }
  if(getdiag) {
    return ( list (
      model = model,
      diag = VARdiag(varxz$y.p, varxz$Z, model$coefficients,
        varxz$ny, varxz$T, (varxz$k+1), varxz$p.max, varxz$dof)
    ))
  } else {
    return (model)
  }
}

.VARXlassocore = function(j, y, x, p, b, lambda, Z, y.spec, x.spec) {
  colIndex = j[1]
  x.to.remove = which(!x.spec[colIndex,])
  y.to.remove = which(!y.spec[colIndex,])
  np = ncol(y) * p
  z.x.to.remove = c()
  z.y.to.remove = c()
  if(length(x.to.remove) > 0) {
    z.x.to.remove = as.vector(sapply(1:b, function(ii) {
      x.to.remove + np + ncol(x)*(ii-1)
    }))
  }
  if(length(y.to.remove) > 0) {
    z.y.to.remove = as.vector(sapply(1:p, function(ii) {
      y.to.remove + ncol(y)*(ii-1)
    }))
  }
  z.to.remove = unlist(c(z.x.to.remove, z.y.to.remove))
  if(length(z.to.remove > 0)) {
    Z.reduced = Z[,-z.to.remove]
  }
  else {
    Z.reduced = Z
  }
       
  if(is.null(lambda)) {
    j.lasso = cv.glmnet(Z.reduced, j[-1])
    j.lasso.s = j.lasso$lambda.min
  } else {
    j.lasso = glmnet(Z.reduced, j[-1])
    j.lasso.s = lambda[colIndex]
  }

  B = rep(0, ncol(Z)+1)  #make vector of coefficients full size
  B.reduced = rep(0, ncol(Z.reduced) + 1)

  j.coef = coef(j.lasso, j.lasso.s)
  B.reduced[(j.coef@i)+1] = j.coef@x
  if(length(z.to.remove) > 0) {
    B[-(z.to.remove+1)] = B.reduced
    return (B)
  }
  else {
    return (B.reduced)
  }
}

VARXlasso = function(y, x, p, b, lambda=NULL, 
  y.spec=matrix(1,nrow=ncol(y),ncol=ncol(y)), 
  x.spec=matrix(1,nrow=ncol(y),ncol=ncol(x)), getdiag=T, numcore=1, ...) {
 
  if(p < 1) stop("p must be a positive integer")
  if(!is.matrix(y)) {
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  varz = VARXZ(y,x,p,b)
  Z = varz$Z
  y.augmented = rbind(1:ncol(y),varz$y.p)
  
  if(numcore==1 | !require(multicore)) {
    var.lasso = apply(y.augmented, 2, .VARXlassocore, y=y,x=x,p=p,b=b,
      lambda=lambda, Z=Z, y.spec=y.spec, x.spec=x.spec)
  } else {
    y.augmented.list = c(unname(as.data.frame(y.augmented)))
    var.lasso = mclapply(y.augmented.list, .VARXlassocore, y=y,x=x,p=p,b=b,
      lambda=lambda, Z=Z, y.spec=y.spec, x.spec=x.spec, mc.cores=numcore, ...)
    var.lasso = matrix(unlist(var.lasso), nrow=(varz$k+1), ncol=ncol(y))
    colnames(var.lasso) = colnames(y)
  }
  rownames(var.lasso) = c('intercept', colnames(Z))
  if(getdiag) {
    return ( list (
      B = var.lasso,
      diag = VARdiag(varz$y.p, Z, var.lasso,
        varz$ny, varz$T, (varz$k+1), varz$p.max, varz$dof)
    ))
  } else {
    return (var.lasso)
  }
}

VARXpredict = function(y, x, p, b, B, xnew, n.ahead=1, threshold) {
  if(p < 1) stop("p must be a positive integer")
  y.pred = matrix(nrow=n.ahead, ncol=ncol(y))
  colnames(y.pred) = colnames(y)
  for(i in 1:n.ahead) {
    Z.ahead.y = as.vector(t(y[((nrow(y)):(nrow(y)-p+1)),]))
    if(b == 0) {
      Z.ahead.x = xnew[1,]
    } else {
      Z.ahead.x = as.vector(t(x[((nrow(x)):(nrow(x)-b+1)),]))
    }
    Z.ahead = c(1, Z.ahead.y, Z.ahead.x)
    y.ahead = Z.ahead %*% B
    if(!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
      y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if(i == n.ahead) break
    y = rbind(y, y.ahead)
    x = rbind(x, xnew[1,])
    xnew = xnew[-1,]
  }
  return (y.pred)
}

VARXZ = function(y, x, p, b) {
  if(p < 1) stop("p must be a positive integer")
  if(is.null(colnames(y))) {
    colnames(y) = sapply(1:ncol(y), function(j) {
      paste('y',j,sep='')
    })
  }
  if(is.null(colnames(x))) {
    colnames(x) = sapply(1:ncol(x), function(j) {
      paste('x',j,sep='')
    })
  }
  
  ny = ncol(y)
  nx = ncol(x)
  T = nrow(y)
  k = ny*p + nx*b
  p.max = max(p,b)
  dof = T - p.max - k

  y.p = y[(1+p.max):T,]
  Z.p = do.call('cbind', (sapply(0:(p-1), function(k) {
    y[(p.max-k):(T-1-k),]
  }, simplify=F)))
  Z.p.names = colnames(Z.p)
  colnames(Z.p) = sapply(1:p, function(i) {
    startIndex = (i-1)*ny + 1
    endIndex = i*ny
    paste(Z.p.names[startIndex:endIndex], '.l', i, sep='')
  })
  if(b == 0) {
    Z.b = x[(1+p):T,]
  } else {
    Z.b = do.call('cbind', (sapply(0:(b-1), function(k) {
      x[(p.max-k):(T-1-k),]
    }, simplify=F)))
    Z.b.names = colnames(Z.b)
    colnames(Z.b) = sapply(1:b, function(i) {
      startIndex = (i-1)*nx + 1
      endIndex = i*nx
      paste(Z.b.names[startIndex:endIndex], '.l', i, sep='')
    })
  }
  Z = cbind(Z.p, Z.b)

  return ( list (
    ny=ny,
    nx=nx,
    T=nrow(y.p),
    k=k,
    p.max=p.max,
    dof=dof,
    y.p=as.matrix(y.p),
    Z=as.matrix(Z)
  ))
}
