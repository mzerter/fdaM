constrain <- function(x,loBnd, upBnd){

  if (x < loBnd)
    {
      x = loBnd
    }
  else if(x > upBnd)
    {
      x = upBnd
    }
  else x
}

