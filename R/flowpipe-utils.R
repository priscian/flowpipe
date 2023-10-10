## Borrowed from package "tidytof" & augmented
## Cf. 'flowCore::arcsinhTransform()': 'asinh(a + b * x) + c'; a = shift about 0, b = scale factor, c = additive constant
#' @export
rev_asinh <- function(
  x,
  shift_factor = 1, # corresponds to parameter 'flowCore::arcsinhTransform(a = 1)'
  scale_factor = 1, # corresponds to parameter 'flowCore::arcsinhTransform(b = 1)'
  additive_constant = 0, # corresponds to parameter 'flowCore::arcsinhTransform(c = 0)',
  a = shift_factor, b = scale_factor, c = additive_constant
)
{
  if (a != shift_factor) shift_factor <- a
  if (b != scale_factor) scale_factor <- b
  if (c != additive_constant) additive_constant <- c

  new_x <- (sinh(x - additive_constant) - shift_factor) / scale_factor

  return(new_x)
}

## usage:
# asinh_test <- function(x, a = 1, b = 1, c = 0) asinh(a + b * x) + c
# asinh_test(1:10, c = 12)
# rev_asinh(asinh_test(1:10, c = 12), c = 12)
