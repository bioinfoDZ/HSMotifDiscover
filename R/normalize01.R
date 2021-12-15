normalize01 <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
