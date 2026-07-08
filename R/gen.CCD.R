#' Generate Central Composite Design (CCD)
#'
#' Generates a Central Composite Design (CCD) according to the specified parameters.
#'
#' @param responses A character vector specifying the name(s) of the response(s).
#' @param factors A character vector specifying the names of the factors.
#' @param ncenter The number of center points to be added in the design.
#' @param randomize Logical, if TRUE, runs are randomized.
#' @param alpha A character among `"rotatable"`, `"spherical"` or `"faced"` that specify the position of the star points.
#'
#' @returns A data.frame where each row corresponds to a run, the first `length(factors)` columns are the factors and the last columns are the responses filled with `NA`.
#'
#'
#'
#' @export
#'
#' @examples
#' gen.CCD(responses=c("R1"),factors=LETTERS[1:4])

gen.CCD=function(responses,factors,ncenter=3,randomize=FALSE,alpha="rotatable"){
  if (!is.character(responses)){
    stop("class(responses) must be cheracter")
  }
  if (is.character(factors)){
    if (length(factors)<2){
      stop("length(factors) must be larger than or equal to 2" )
    }
  }else{
    stop("class(factors) must be cheracter")
  }
  if (is.numeric(ncenter) | is.integer(ncenter)){
    if (ncenter<0 | ncenter%%1!=0){
      stop("ncenter must be a positive integer")
    }
  }else{
    stop("class(ncenter) must be numeric or integer")
  }
  if (!is.logical(randomize)){
    stop("class(randomize) must be logical")
  }
  if (is.character(alpha)){
    if (!alpha%in%c("rotatable","spherical","faced")){
      stop("alpha must be rotatable, spherical or faced")
    }
  }else{
    stop("class(alpha) must be character")
  }
  D=expand.grid(replicate(length(factors),c(-1,1),simplify = FALSE)) ; colnames(D)=factors
  if (alpha=="rotatable"){
    alpha=nrow(D)^(1/4)
  }else if (alpha=="spherical"){
    alpha=sqrt(length(factors))
  }else{
    alpha=1
  }
  for (j in 1:ncol(D)){
    ajout=matrix(0,2,ncol(D))
    ajout[,j]=c(-alpha,alpha) ; colnames(ajout)=colnames(D)
    D=rbind(D,ajout)
  }
  if (ncenter>0){
    zero=as.data.frame(matrix(0,ncenter,ncol(D))) ; colnames(zero)=colnames(D)
    D=rbind(D,zero)
  }
  if (randomize){
    D=D[sample(1:nrow(D)),]
  }
  re=cbind.data.frame(D,matrix(NA,nrow(D),length(responses))) ; colnames(re)[-c(1:ncol(D))]=responses
  class(re)=c("DoETeachR","data.frame")
  attr(re,"colonne")=c(rep("factor",length(factors)),rep("response",length(responses)))
  attr(re,"type")="CCD"
  return(re)
}
