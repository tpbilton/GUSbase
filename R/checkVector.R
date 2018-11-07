

checkVector <- function(x, type="pos_integer", minv=0, maxv=NULL, equal=TRUE){
  if(type=="pos_integer" & !is.null(maxv))
    return(!is.vector(x) || !is.numeric(x) || any(is.na(x)) || !all(x == round(x)) || any(x < minv) || any(x > maxv))
  else if(type=="pos_integer" & is.null(maxv))
    return(!is.vector(x) || !is.numeric(x) || any(is.na(x)) || !all(x == round(x)) || any(x < minv) || is.infinite(x))
  else if(type=="pos_numeric" & !is.null(maxv)){
    if(equal)
      return(!is.vector(x) || !is.numeric(x) || any(is.na(x)) || any(x <= minv) || any(x >= maxv))
    else
      return(!is.vector(x) || !is.numeric(x) || any(is.na(x)) || any(x < minv) || any(x > maxv))
  }
  else if(type=="pos_numeric" & is.null(maxv)){
    if(equal)
      return(!is.vector(x) || !is.numeric(x) || any(is.na(x)) || any(x <= minv) || is.infinite(x))
    else
      return(!is.vector(x) || !is.numeric(x) || any(is.na(x)) || any(x < minv) || is.infinite(x))
  } else if(type == "one_logical")
    return(!is.logical(x) || any(is.na(x)) || length(x) != 1)
}
