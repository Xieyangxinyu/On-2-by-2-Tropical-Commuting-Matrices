## ---------------------------
##
## Script name: Extreme vertices of 2 by 2 tropical commuting matrices
##
## Purpose of script: This script aims to give an easy implementation for the 
##                      visualisation of the Commuting Polyhedral Cone in 
##                      Barycentric coordinates 
##
## Author: Yangxinyu Xie
##
## Date Created: 2019-June-09
##
## Copyright (c) Yangxinyu Xie, 2019
## Email: yx4247@utexas.edu
##
## ---------------------------
## Reference: Xie, Y(2019). On 2 by 2 Tropical Commuting Matrices.
## ---------------------------

# Make sure you install the following packages if you haven't done so
#install.packages("geometry")
#install.packages("retistruct")

# Import Libraries
library(retistruct)
library(geometry)

# Input: 
## beta: a list of entries of the extreme vertices of the 
##       tropical polyhedral cone defined by the two-sided
##       system with commuting constraints
## vertices_label: the labels of the coordinates
# Output: Give a Barycentric visualisation if a_11 != a_22
draw_Barycentric <- function(beta, vertices_label){
  
  ## Define simplex in 2D
  X <- rbind(c(0, 0),
             c(1, 0),
             c(1/2, sqrt(3)/2))
  ## Draw the simplex
  trimesh(rbind(1:3), X)
  ## Label coordinates
  text(X[,1], X[,2] + 0.02, vertices_label) # Label vertices
  ## Calculate the Cartesian coordinates
  P <- bary2cart(X, beta)
  ## Attach the intersection point to P
  P <- rbind(P, line.line.intersection(c(1, 0), c(P[2,1], P[2,2]), c(0,0), c(P[3,1], P[3,2]), interior.only = FALSE))
  coordinate <- bary2cart(X, rbind(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  ))
  ## Draw disect the area into subareas
  for (i in 2:4){
    segments(x0 = P[i, 1], y0 = P[i, 2], x1 = coordinate[1,1], y1 = coordinate[1,2], lty = 2, lwd = 0.5)
    segments(x0 = P[i, 1], y0 = P[i, 2], x1 = coordinate[2,1], y1 = coordinate[2,2], lty = 2, lwd = 0.5)
    segments(x0 = P[i, 1], y0 = P[i, 2], x1 = coordinate[3,1], y1 = coordinate[3,2], lty = 2, lwd = 0.5)
  }
  ## Shade the cell within the Tropical Polyhedral Cone
  polygon(c(P[c(1,2,5,3), 1]),c(P[c(1,2,5,3), 2]),col="grey")
  ## Circle the extreme points
  points(P)
  
  ## Draw the edges of the Tropical Polyhedral Cone
  for (i in 2:4){
    segments(x0 = P[i, 1], y0 = P[i, 2], x1 = P[5,1], y1 = P[5,2], lty = 1, lwd = 2)
  }
  segments(x0 = P[1, 1], y0 = P[1, 2], x1 =  P[2,1], y1 = P[2,2], lty = 1, lwd = 2)
  segments(x0 = P[1, 1], y0 = P[1, 2], x1 = P[3,1], y1 = P[3,2], lty = 1, lwd = 2)
}

# Input: 
## beta: a two by two matrix A
# Output: a list of matrices whose entries are the extreme
##        vertices
get_B <- function(A){
  ## Error Handling
  if (nrow(A) != 2 || ncol(A) != 2){
    stop("Only 2 by 2 matrices are supported")
  }
  
  ## Get the entries of A
  a_11 = A[1,1]
  a_12 = A[1,2]
  a_21 = A[2,1]
  a_22 = A[2,2]
  if (a_11 == -Inf || a_11 == Inf
      || a_12 == -Inf || a_12 == Inf
      || a_21 == -Inf || a_21 == Inf
      || a_22 == -Inf || a_22 == Inf){
    stop("Matrix A cannot have Infinite values")
  }
  
  a_star = min(a_12 - a_11, a_22 - a_21, a_11 - a_21, a_12-a_22)
  a_prime = a_star + a_21 - a_12
  
  ## All potential betas as described in the reference
  beta_pool <- list(matrix(data = c(0, -Inf, -Inf, 0), nrow = 2),
                    matrix(data = c(0, a_star, -Inf, 0), nrow = 2),
                    matrix(data = c(0, -Inf, a_prime, 0), nrow = 2),
                    matrix(data = c(0, a_12 - a_11, a_21 - a_11, -Inf), nrow = 2),
                    matrix(data = c(-Inf, 0, a_21 - a_12, a_22-a_12), nrow = 2),
                    matrix(data = c(-Inf, 0, a_21 - a_12, 0), nrow = 2))
  
  if (a_11 > a_22){
    beta = beta_pool[c(1,2,3,4)]
    ## Get the Barycentric coordinates
    beta_Bary <- rbind(c(0, 0, 1),
                  c(exp(a_star)/(exp(a_star) + 1), 0, 1/(exp(a_star) + 1)),
                  c(0, exp(a_prime)/(exp(a_prime) + 1), 1/(exp(a_prime) + 1)),
                  c(exp(a_12 - a_11)/(exp(a_12 - a_11) + exp(a_21 - a_11)), exp(a_21 - a_11)/(exp(a_12 - a_11) + exp(a_21 - a_11)), 0)
    )
    draw_Barycentric(beta_Bary, c('b12', 'b21', 'b22'))
  }else if (a_11 < a_22){
    beta = beta_pool[c(1,2,3,5)]
    beta_Bary <- rbind(c(0, 0, 1),
                       c(exp(a_star)/(exp(a_star) + 1), 0, 1/(exp(a_star) + 1)),
                       c(0, exp(a_prime)/(exp(a_prime) + 1), 1/(exp(a_prime) + 1)),
                       c(exp(a_12 - a_22)/(exp(a_12 - a_22) + exp(a_21 - a_22)), exp(a_21 - a_22)/(exp(a_12 - a_22) + exp(a_21 - a_22)), 0)
    )
    draw_Barycentric(beta_Bary, c('b12', 'b21', 'b11'))
  }
  else{
    beta = beta_pool
    print("Visualisation for the case when a_11 = a_22 is not supported")
  }
  return(beta)
}

#Example
A <- matrix(data = c(1.2, 1.3, 1.2, 1.1), nrow = 2)
print(get_B(A))

A <- matrix(data = c(-1, 1.3, 1.2, 1.1), nrow = 2)
print(get_B(A))

