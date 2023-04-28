## Lifted from package "tidytof"
rev_asinh <- function(
  x,
  shift_factor,
  scale_factor
)
{
  new_x <- (sinh(x) - shift_factor)/scale_factor

  return(new_x)
}
