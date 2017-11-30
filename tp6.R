### TP6 ###
##Charles-Alexandre Bernard
## Nicolas H?bert
## Marc-Andr? Noel


#1


x<-scan("travail_pratique_6.dat")


nr.mult <- function(FUN, gradiant, start, TOL = 1E-6,
                    maxiter = 100, echo = FALSE, ...)
{
    if (echo)
        expr <- expression(print(start <- start - adjust))
    else
        expr <- expression(start <- start - adjust)

    i <- 0                           # initialisation du compteur

    repeat
    {
        adjust <- solve(gradiant(start, ...), FUN(start, ...))

        if (max(abs(adjust)) < TOL)
            break

        if (maxiter < (i <- i + 1))
            stop("Maximum number of iterations reached
                  without convergence")

        eval(expr)
    }
    list(roots = start - adjust, nb.iter = i)
}

f <- function(p,x){
  a <- p[1]
  l <- p[2]
  c(n/a+n*log(l)-sum(log(x+l)),n*a/l-(a+1)*sum(1/(x+l)))
}

fp <- function(p,x){
  a <- p[1]
  l <- p[2]
  cbind(c(-n/a^2,n/l-sum( 1/(x+l))),c(n/l-sum(1/(x+l)),-n*a/l^2+(a+1)*sum((x+l)^-2)))
}

(r <- nr.mult(f,fp,start=c(2,1000),TOL = 1E-6,maxiter = 10000, echo = FALSE, x=x))


#2
#a)
a <- matrix(c(-8.64,0.101,1.997,-1.095),nrow=2)
(vp <- eigen(a))
p <- vp$vectors
D <- diag(exp(vp$values))
(p1 <- solve(p))

k<- p%*%D%*%p1

pie<- c(0.5614, 0.4386)

j <- pie%*%k
e <- matrix(c(1,1))
1- (j %*% e)

#b)
pphtype(1,pie,a)

salut
