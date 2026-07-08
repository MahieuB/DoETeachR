#' Augment a design of experiments from \pkg{DoETeachR}
#'
#' Augments a design of experiments from \pkg{DoETeachR} depending on its type. See details.
#'
#' @param design A design from \code{\link[DoETeachR]{gen.FrFD}}, \code{\link[DoETeachR]{gen.FFD}} or \code{\link[DoETeachR]{gen.CCD}}.
#' @param randomize Logical, if TRUE, added runs are randomized.
#' @param star.res A positive integer giving the resolution at which \code{\link[DoETeachR]{AugmentDesign}} adds star points rather than a fraction for Fractional Factorial Designs.
#' @param alpha A character among `"rotatable"`, `"spherical"` or `"faced"` that specify the position of the star points.
#'
#'
#' @returns A data.frame where each row corresponds to a run, the first `length(factors)` columns are the factors and the last columns are the responses filled with `NA` for the new runs.
#'
#' @details For FFD it adds star points. For FrFD it adds a fraction up to a \bold{star.res} resolution design and then adds star points. For CCD, it replicates the design.
#'
#' @import stats
#'
#' @export
#'
#' @examples
#' design=gen.FFD(responses=c("R1"),factors=LETTERS[1:4])
#' AugmentDesign(design)
#'
#' design=gen.FrFD(responses=c("R1"),factors=LETTERS[1:4])
#' AugmentDesign(design)
#'
#' design=gen.CCD(responses=c("R1"),factors=LETTERS[1:4])
#' AugmentDesign(design)

AugmentDesign=function(design,randomize=FALSE,star.res=5,alpha="rotatable"){
  if (!inherits(design,"DoETeachR")){
    stop("design must be an ExpeRimenteR object")
  }
  if (!is.logical(randomize)){
    stop("class(randomize) must be logical")
  }
  if (is.integer(star.res) | is.numeric(star.res)){
    if (star.res<3){
      stop("star.res must be larger than or equal to 3")
    }
  }else{
    stop("class(star.res) must be numeric or integer")
  }
  if (is.character(alpha)){
    if (!alpha%in%c("rotatable","spherical","faced")){
      stop("alpha must be rotatable, spherical or faced")
    }
  }else{
    stop("class(alpha) must be character")
  }
  if (attr(design,"type")=="FrFD"){
    nfac=sum(attr(design,"colonne")=="factor") ; nrepon=sum(attr(design,"colonne")=="response")
    lesi=attr(design,"alias") ; minr=min(sapply(strsplit(lesi,"[:]"),length))
    if (minr<star.res){
      if (length(which(rowSums(design[,1:nfac]==0)==nfac))>0){
        design.sans0=design[-which(rowSums(design[,1:nfac]==0)==nfac),]
        nsans0=nrow(design.sans0)
      }else{
        nsans0=nrow(design)
      }
      lesimoins=lesi[sapply(strsplit(lesi,"[:]"),length)==minr] ; lesiplus=lesi[!sapply(strsplit(lesi,"[:]"),length)==minr]
      D.full=expand.grid(replicate(nfac,c(-1,1),simplify = FALSE)) ; colnames(D.full)=colnames(design[1:nfac])
      mod.mat.full=model.matrix(as.formula(paste("~",paste(rep(".",nfac),collapse = "*"),sep="")),D.full)
      if (length(lesiplus)>0 & length(lesimoins)>0){
        ajout=D.full[apply(mod.mat.full[,lesimoins,drop=FALSE]==-1,1,all) & apply(mod.mat.full[,lesiplus,drop=FALSE]==1,1,all),]
      }else if (length(lesiplus)==0 & length(lesimoins)>0){
        ajout=D.full[apply(mod.mat.full[,lesimoins,drop=FALSE]==-1,1,all),]
      }
      droite=matrix(NA,nrow(ajout),nrepon)
      ajout=cbind(ajout,droite) ; colnames(ajout)=colnames(design)
      if (randomize){
        ajout=ajout[sample(1:nrow(ajout)),]
      }
      re=rbind(design,ajout)
      class(re)=c("DoETeachR","data.frame")
      attr(re,"colonne")=attr(design,"colonne")
      if (length(lesimoins)==1 & length(lesiplus)==0){
        vec.alias=NULL
      }else{
        if (length(lesimoins)>1){
          vec.alias=NULL
          for (j in 2:length(lesimoins)){
            product=paste(lesimoins[1],lesimoins[j],sep=":") ; split.product=strsplit(product,"[:]")[[1]]
            stay=names(table(split.product))[table(split.product)==1]
            vec.alias=c(vec.alias,paste(stay,collapse = ":"))
          }
          vec.alias=c(vec.alias,lesiplus)
        }else{
          vec.alias=lesiplus
        }
      }
      attr(re,"alias")=vec.alias
      if ((nrow(ajout)+nsans0)==2^nfac){
        attr(re,"type")="FFD"
      }else{
        attr(re,"type")="FrFD"
      }
    }else{
      ajout=NULL
      if (alpha=="rotatable"){
        if (length(which(rowSums(design[,1:nfac]==0)==nfac))>0){
          design.sans0=design[-which(rowSums(design[,1:nfac]==0)==nfac),]
          nsans0=nrow(design.sans0)
        }else{
          nsans0=nrow(design)
        }
        alpha=nsans0^(1/4)
      }else if (alpha=="spherical"){
        alpha=sqrt(nfac)
      }else{
        alpha=1
      }
      for (j in 1:nfac){
        ajout.loc=matrix(0,2,nfac)
        ajout.loc[,j]=c(-alpha,alpha)
        ajout=rbind(ajout,ajout.loc)
      }
      droite=matrix(NA,nrow(ajout),nrepon)
      ajout=cbind(ajout,droite) ; colnames(ajout)=colnames(design)
      if (randomize){
        ajout=ajout[sample(1:nrow(ajout)),]
      }
      re=rbind(design,ajout)
      class(re)=c("DoETeachR","data.frame")
      attr(re,"colonne")=attr(re,"colonne")
      attr(re,"type")="CCD"
    }
  }else if (attr(design,"type")=="FFD"){
    nfac=sum(attr(design,"colonne")=="factor") ; nrepon=sum(attr(design,"colonne")=="response")
    ajout=NULL
    if (alpha=="rotatable"){
      if (length(which(rowSums(design[,1:nfac]==0)==nfac))>0){
        design.sans0=design[-which(rowSums(design[,1:nfac]==0)==nfac),]
        nsans0=nrow(design.sans0)
      }else{
        nsans0=nrow(design)
      }
      alpha=nsans0^(1/4)
    }else if (alpha=="spherical"){
      alpha=sqrt(nfac)
    }else{
      alpha=1
    }
    for (j in 1:nfac){
      ajout.loc=matrix(0,2,nfac)
      ajout.loc[,j]=c(-alpha,alpha)
      ajout=rbind(ajout,ajout.loc)
    }
    droite=matrix(NA,nrow(ajout),nrepon)
    ajout=cbind(ajout,droite) ; colnames(ajout)=colnames(design)
    if (randomize){
      ajout=ajout[sample(1:nrow(ajout)),]
    }
    re=rbind(design,ajout)
    class(re)=c("DoETeachR","data.frame")
    attr(re,"colonne")=attr(re,"colonne")
    attr(re,"type")="CCD"
  }else{
    nfac=sum(attr(design,"colonne")=="factor") ; nrepon=sum(attr(design,"colonne")=="response")
    ajout=design[1:nfac]
    droite=matrix(NA,nrow(ajout),nrepon)
    ajout=cbind(ajout,droite) ; colnames(ajout)=colnames(design)
    if (randomize){
      ajout=ajout[sample(1:nrow(ajout)),]
    }
    re=rbind(design,ajout)
    class(re)=c("DoETeachR","data.frame")
    attr(re,"colonne")=attr(re,"colonne")
    attr(re,"type")=attr(re,"type")
  }
  rownames(re)=as.character(1:nrow(re))
  return(re)
}
