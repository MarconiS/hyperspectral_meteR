# utility functions
#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
pls_transform <- function(Xt, Yt, nComp){
  library(MASS)
  X_ave = apply(Xt, 2, mean)
  X_sd = apply(Xt, 2, sd)
  Y_ave = apply(Yt[-1], 2, mean)
  Y_sd = apply(Yt[-1], 2, sd)
  X <- scale(Xt, center = X_ave, scale = X_sd)
  Y <- scale(Yt[-1], center = Y_ave, scale = Y_sd)
  C = t(X) %*% Y
  C_uv = svd(C)
  C_uv = svd_flip(U = C_uv$u, V = C_uv$v)
  V = t(C_uv$v)
  x_scores = X %*% C_uv$u
  y_scores = Y %*% V
  x_weights = C_uv$u
  y_weights = V
  wwgts <- list(x_weights = x_weights, y_weights = y_weights,
                x_scores = x_scores, y_scores = y_scores, 
                x_center = X_ave, x_scale = X_sd,
                y_center = Y_ave, y_scale = Y_sd)
  return(wwgts)
}

# utility functions
#' convert utmZone into corresponding EPSG
#'
#'
#' @inheritParams str_detect
#' @return A list of dataframe
#' @import rotl ape
#' @examples
#' @importFrom magrittr "%>%"
svd_flip <- function(U, V){
  #  max_abs_cols = np.argmax(np.abs(u), axis=0)
  #Adjusts the columns of u and the rows of v such that the loadings in the
  #columns in u that are largest in absolute value are always positive.
  max_abs_cols = unlist(apply(U, 2, function(x) {which(abs(x) == max(abs(x)), arr.ind = TRUE)})) #%>%
  signs = unlist(lapply(1:dim(U)[2], function(x) (sign(U[max_abs_cols[x], x]))))
  u = U * signs
  v = V * signs
  return(list(u = u, v = v))
}

# 4) rotations from input space to transformed space (scores)
# T = X W(P'W)^-1 = XW* (W* : p x k matrix)
# U = Y C(Q'C)^-1 = YC* (W* : q x k matrix)

