VAR.Z = function(y, p) {
  if(p < 1) stop("p must be a positive integer")
      
  if(is.null(colnames(y))) {
    colnames(y) = sapply(1:ncol(y), function(j) {
      paste('y',j,sep='')
    })
  }

  n = ncol(y)
  T = nrow(y)
  k = n*p
  dof = T - p - k
        
  y.p = y[(1+p):T,]
  Z = do.call('cbind', (sapply(0:(p-1), function(k) {
    y[(p-k):(T-1-k),]
  }, simplify=F)))
  Z.names = colnames(Z)
  colnames(Z) = sapply(1:p, function(i) {
    startIndex = (i-1)*n + 1
    endIndex = i*n
    paste(Z.names[startIndex:endIndex], '.l', i, sep='')
  })

  return ( structure ( list (
    n=n,
    T=nrow(y.p),
    k=k,
    dof=dof,
    y.p=as.matrix(y.p),
    Z=as.matrix(Z),
    y.orig=y,
    p = p
  ), class="fastVAR.VARZ"))
}



VARX.Z = function(y, x, p, b) {
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
    b = b
  ), class="fastVAR.VARXZ"))
}
