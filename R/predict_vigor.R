predict_vigor <- function(object, newX){

  #newX is a list of X to be predicted
  #The length and order of newX are the same as those of ETA used for training
  #For BLUP, X is the n1 x n2 matrix of relationship matrix
  ##where n1 and n2 are the numbers of samples in test and training data

  #Transform formula if included
  Class <- lapply(newX, class)
  if(any(Class == "formula")){
    if(!any(names(newX) == "data"))
      stop("Provide data when fixed effects are specified with formula")

    Data <- newX$data
    w <- which(Class == "formula")
    d <- which(names(newX) == "data")
    newX[[w]] <- model.matrix(newX[[w]], data = Data)
    newX <- newX[-d]
  }

  #Number of samples in newX
  N <- unlist(lapply(newX, nrow))
  if(length(unique(N))>1)
    stop("Number of samples is inconsistent between methods in newX")

  N <- N[1]

  if(object$AddIntercept){
    if(!is.list(newX) | length(newX)+1 != length(object$ETA))
      stop("newX should be a list with the same length as ETA used for training")

    newX<-c(newX, list(matrix(1, N, ncol=1)))#Intercept

  }else{
    if(!is.list(newX) | length(newX) != length(object$ETA))
      stop("newX should be a list with the same length as ETA used for training")
  }

  Nm <- length(object$ETA)
  Prediction <- as.list(numeric(Nm))

  for(method in 1:Nm){
    if(is.null(newX[[method]])){
      Prediction[[method]] <- NULL
    }else{
      if(!is.null(object$ETA[[method]]$Beta)){
        #Other than BLUP
        if(ncol(newX[[method]]) != length(object$ETA[[method]]$Beta))
          stop("Number of covariates in method", method,
               "is incosistent with object")
        Prediction[[method]] <-
          newX[[method]] %*%
          matrix(object$ETA[[method]]$Beta,ncol=1)
      }else{
        #BLUP
        if(ncol(newX[[method]]) != length(object$ETA[[method]]$U))
          stop("Number of levels in method", method,
               " is incosistent with object")
        Prediction[[method]] <-
          newX[[method]] %*%
          object$ETA[[method]]$iK %*%
          matrix(object$ETA[[method]]$U, ncol=1)
      }
    }
  }

  Prediction.sum <- numeric(N)
  for(method in 1:Nm){
    if(!is.null(newX[[method]]))
      Prediction.sum <- Prediction.sum + as.vector(Prediction[[method]])
  }

  Prediction.sum
}
