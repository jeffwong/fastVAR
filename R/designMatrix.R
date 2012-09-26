VAR.Z = function(y, p, intercept=F) {
  if(p < 1) stop("p must be a positive integer")
      
  if(is.null(colnames(y))) {
    colnames(y) = sapply(1:ncol(y), function(j) {
      paste('y',j,sep='')
    })
  }

  n = ncol(y) #numcols in response
  T = nrow(y) #numrows in response
  k = n*p     #numcols in design matrix
  if(intercept) k = k+1
  dof = T - p - k
        
  y.p = y[(1+p):T,]
  Z = do.call('cbind', (sapply(0:(p-1), function(i) {
    y[(p-i):(T-1-i),]
  }, simplify=F)))
  Z.names = colnames(Z)
  colnames(Z) = as.vector(sapply(1:p, function(i) {
    startIndex = (i-1)*n + 1
    endIndex = i*n
    paste(Z.names[startIndex:endIndex], '.l', i, sep='')
  }))

  if(intercept) {
    Z = cbind(1, Z)
    colnames(Z)[1] = "intercept"
  }

  return ( structure ( list (
    n=n,
    T=nrow(y.p),
    k=k,
    dof=dof,
    y.p=as.matrix(y.p),
    Z=as.matrix(Z),
    y.orig=y,
    p = p,
    intercept = intercept
  ), class="fastVAR.VARZ"))
}



VARX.Z = function(y, x, p, b, intercept=F) {
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
  if(intercept) k = k+1
  p.max = max(p,b)
  dof = T - p.max - k

  y.p = y[(1+p.max):T,]
  Z.p = do.call('cbind', (sapply(0:(p-1), function(i) {
    y[(p.max-i):(T-1-i),]
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
  if(intercept) {
      Z = cbind(1, Z)
      colnames(Z)[1] = "intercept"
  }

  return ( structure ( list (
    ny=ny,
    nx=nx,
    T=nrow(y.p),
    k=k,
    p.max=p.max,
    dof=dof,
    y.p=as.matrix(y.p),
    Z=as.matrix(Z),
    y.orig = y,
    x.orig = x,
    p = p,
    b = b,
    intercept = intercept
  ), class="fastVAR.VARXZ"))
}
