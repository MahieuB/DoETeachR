#' Generate Fractional Factorial Design (FrFD)
#'
#' Generates a Fractional Factorial Design (FrFD) according to the specified parameters.
#'
#' @param responses A character vector specifying the name(s) of the response(s).
#' @param factors A character vector specifying the names of the factors.
#' @param frac An integer power of 2. The fraction of the full factorial design. For example `frac=2` corresponds to half the full factorial design i.e. `2^(length(factors)-1)` runs.
#' @param ncenter The number of center points to be added in the design.
#' @param randomize Logical, if TRUE, runs are randomized.
#'
#' @returns A data.frame where each row corresponds to a run, the first `length(factors)` columns are the factors and the last columns are the responses filled with `NA`.
#'
#' @import stats
#'
#' @export
#'
#' @examples
#' gen.FrFD(responses=c("R1"),factors=LETTERS[1:4])

gen.FrFD=function(responses,factors,frac=2^1,ncenter=0,randomize=FALSE){
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
  if (is.numeric(frac) | is.integer(frac)){
    if (log2(frac)<0 | log2(frac)%%1!=0){
      stop("frac must be a power of 2")
    }
  }else{
    stop("class(frac) must be numeric or integer")
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
  D.gen=expand.grid(replicate(length(factors)-log2(frac),c(-1,1),simplify = FALSE)) ; colnames(D.gen)=factors[1:ncol(D.gen)]
  mod.mat=model.matrix(as.formula(paste("~",paste(rep(".",ncol(D.gen)),collapse = "*"),sep="")),D.gen)
  vec.alias=NULL ; vec.candidat=1:ncol(mod.mat)
  for (f in 1:log2(frac)){
    D.gen=cbind(D.gen,mod.mat[,vec.candidat[length(vec.candidat)]]) ; colnames(D.gen)[ncol(D.gen)]=factors[length(factors)-f+1]
    vec.alias=c(vec.alias,paste(colnames(mod.mat)[vec.candidat[length(vec.candidat)]],factors[length(factors)-f+1],sep=":"))
    vec.candidat=vec.candidat[-length(vec.candidat)]
  }
  if (length(vec.alias)>1){
    minr=min(sapply(strsplit(vec.alias,"[:]"),length))
    taille=minr
    for (i.1 in 1:(length(vec.alias)-1)){
      for (i.2 in (i.1+1):length(vec.alias)){
        prodi=paste(vec.alias[i.1],vec.alias[i.2],sep=":") ; split.prodi=strsplit(prodi,"[:]")[[1]]
        taille=min(taille,length(names(table(split.prodi))[table(split.prodi)==1]))
      }
    }
    decal=1
    while (taille<minr){
      D.gen=expand.grid(replicate(length(factors)-log2(frac),c(-1,1),simplify = FALSE)) ; colnames(D.gen)=factors[1:ncol(D.gen)]
      mod.mat=model.matrix(as.formula(paste("~",paste(rep(".",ncol(D.gen)),collapse = "*"),sep="")),D.gen)
      vec.alias=NULL ; vec.candidat=1:(ncol(mod.mat)-decal)
      for (f in 1:log2(frac)){
        D.gen=cbind(D.gen,mod.mat[,vec.candidat[length(vec.candidat)]]) ; colnames(D.gen)[ncol(D.gen)]=factors[length(factors)-f+1]
        vec.alias=c(vec.alias,paste(colnames(mod.mat)[vec.candidat[length(vec.candidat)]],factors[length(factors)-f+1],sep=":"))
        vec.candidat=vec.candidat[-length(vec.candidat)]
      }
      minr=min(sapply(strsplit(vec.alias,"[:]"),length))
      taille=minr
      for (i.1 in 1:(length(vec.alias)-1)){
        for (i.2 in (i.1+1):length(vec.alias)){
          prodi=paste(vec.alias[i.1],vec.alias[i.2],sep=":") ; split.prodi=strsplit(prodi,"[:]")[[1]]
          taille=min(taille,length(names(table(split.prodi))[table(split.prodi)==1]))
        }
      }
      decal=decal+1
    }
  }
  D=D.gen[,factors[1:length(factors)]] ; colnames(D)=factors
  if (ncenter>0){
    zero=as.data.frame(matrix(0,ncenter,ncol(D))) ; colnames(zero)=colnames(D)
    D=rbind(D,zero)
  }
  if (randomize){
    D=D[sample(1:nrow(D)),]
  }
  re=cbind.data.frame(D,matrix(NA,nrow(D),length(responses))) ; colnames(re)[-c(1:ncol(D))]=responses
  rownames(re)=as.character(1:nrow(re))
  class(re)=c("DoETeachR","data.frame")
  attr(re,"colonne")=c(rep("factor",length(factors)),rep("response",length(responses)))
  attr(re,"type")="FrFD"
  attr(re,"alias")=rev(vec.alias)
  return(re)
}
