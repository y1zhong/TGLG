# function for adjust the acceptance rate
adjust_acceptance=function(accept,sgm,target = 0.5){
  y = 1. + 1000.*(accept-target)*(accept-target)*(accept-target)
  if (y < .9)
    y = .9
  if (y > 1.1)
    y = 1.1
  sgm = sgm* y
  return(sgm)
}

# approximate the threhold function
appro2 = function(x, epison = 10^-8) {
  temp = 0.5*(1+2/pi*atan(x/epison))
  return(temp)
}

# appraoximate the derivative of the threshold function
dappro2 = function(x, epison = 10^-8) {
  temp = 2/(pi*(epison+x^2/epison))
  return(temp)
}
