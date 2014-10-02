fdParcheck <- function(fdParobj){
  if (is.fdPar(fdParobj)) {
    
    if (is.fd(fdParobj) || is.basis(fdParobj)) {
      
      fdParobj = fdPar(fdParobj)
    }
    else{
      stop('FDPAROBJ is not a functional parameter object, not a functional data object, and not a basis object.')
    }
      
  }
}