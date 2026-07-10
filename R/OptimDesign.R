#' Compute Optimal Factor Values for a Response Target of \pkg{DoETeachR} Objects
#'
#' Computes the optimal values of the factors to reach \bold{target.optim} of the response given \bold{formula}.
#'
#' @param formula A formula that defines the model. Only order two effects are supported at maximum.
#' @param design A design from \code{\link[DoETeachR]{gen.FrFD}}, \code{\link[DoETeachR]{gen.FFD}} or \code{\link[DoETeachR]{gen.CCD}}.
#' @param bounds Either NULL (Default) or a numeric vector of length two defining the range of values that must be considered for the factors. If NULL, the range of values for the factors are the minimum and maximum values of the plotted factors present in \bold{formula}.
#' @param target.optim Either `"min"`, `"max"` or a numeric value corresponding to a target value of the response indicating whether the objective is to minimize, maximize or maintain the value of the response at the target respectively.
#'
#' @returns A list if length two where the first element is the optimal values of the factors and the second element is the value of the response obtained at the optimal values of the factors.
#'
#' @import stats
#'
#' @export
#'
#' @examples
#' design=gen.CCD(responses=c("R1"),factors=LETTERS[1:4])
#' design$R1=runif(nrow(design))
#' OptimDesign(R1~.*.+I(A^2)+I(B^2)+I(C^2)+I(D^2),design)

OptimDesign=function(formula,design,bounds=NULL,target.optim="min"){
  if (inherits(formula,"formula")){
    droite=as.character(formula)[3]
    if (any(gregexpr("poly\\(",droite)[[1]]>0)){
      stop("formula must no contain poly(")
    }
    if (length(gregexpr("[.]",droite)[[1]])>2){
      stop("formula must no contain effect of order larger than 2 and more than two dots (.)")
    }
    if (any(gregexpr("\\^3",droite)[[1]]>0) | any(gregexpr("\\^4",droite)[[1]]>0) | any(gregexpr("\\^5",droite)[[1]]>0)){
      stop("formula must no contain effect of order larger than 2")
    }
  }else{
    stop("class(formula) must be formula")
  }
  if (!inherits(design,"DoETeachR")){
    stop("class(design) must be DoETeachR")
  }
  mod=lm(formula,design) ; x=model.matrix(mod) ; b=coef(mod)
  verif.alias=alias(mod)
  if (!is.null(verif.alias$Complete)){
    stop("There are aliased term in the model")
  }
  if (any(sapply(strsplit(colnames(x),"[:]"),length)>2)){
    stop("formula must no contain effect of order larger than 2")
  }
  factor.name=colnames(design)[attr(design,"colonne")=="factor"] ; factor.name=factor.name[apply(sapply(factor.name, gregexpr,text=colnames(x))>0,2,any)]
  if (is.numeric(bounds) | is.null(bounds)){
    if (!is.null(bounds)){
      if (length(bounds)!=2){
        stop("length(bounds) must be 2 when class(bounds) is numeric")
      }
    }
  }else{
    stop("class(bounds) must be numeric or NULL")
  }
  if (is.character(target.optim) | is.numeric(target.optim)){
    if (is.character(target.optim)){
      if (!target.optim%in%c("min","max")){
        stop("target.optim must be min or max when class(target.optim) is character")
      }
    }
  }else{
    stop("class(target.optim) must be character or numeric")
  }
  if (is.null(bounds)){
    bounds=c(min(x),-min(x))
  }
  if (is.character(target.optim)){
    fopt=function(vec){
      pseudo.D=as.data.frame(t(as.matrix(vec))) ; colnames(pseudo.D)=factor.name ;pseudo.formula=as.formula(paste(as.character(formula)[1],as.character(formula)[3],sep=""))
      pseudo.mat=model.matrix(pseudo.formula,pseudo.D)
      y=as.numeric(as.matrix(pseudo.mat)%*%b[colnames(pseudo.mat)])
      return(y)
    }
    if (target.optim=="max"){
      y.opti=-Inf
      x.pure=unique(x[(x%*%b)>=quantile(x%*%b,probs = 0.60,type=5),factor.name,drop=FALSE])
      for (i in 1:nrow(x.pure)){
        opti.i=optim(x.pure[i,],fopt,method = "L-BFGS-B",lower=rep(bounds[1],length(factor.name)),upper=rep(bounds[2],length(factor.name)),control = list(fnscale=-1))
        if (opti.i$value>y.opti){
          opti=opti.i
          y.opti=opti$value
        }
      }
    }else{
      y.opti=Inf
      x.pure=unique(x[(x%*%b)<=quantile(x%*%b,probs = 0.40,type=5),factor.name,drop=FALSE])
      for (i in 1:nrow(x.pure)){
        opti.i=optim(x.pure[i,],fopt,method = "L-BFGS-B",lower=rep(bounds[1],length(factor.name)),upper=rep(bounds[2],length(factor.name)),control = list(fnscale=1))
        if (opti.i$value<y.opti){
          opti=opti.i
          y.opti=opti$value
        }
      }
    }
  }else{
    fopt=function(vec){
      pseudo.D=as.data.frame(t(as.matrix(vec))) ; colnames(pseudo.D)=factor.name ;pseudo.formula=as.formula(paste(as.character(formula)[1],as.character(formula)[3],sep=""))
      pseudo.mat=model.matrix(pseudo.formula,pseudo.D)
      y=as.numeric(as.matrix(pseudo.mat)%*%b[colnames(pseudo.mat)])
      return(sqrt((y-target.optim)^2))
    }
    opti=optim(rep(0,length(factor.name)),fopt,method = "L-BFGS-B",lower=rep(bounds[1],length(factor.name)),upper=rep(bounds[2],length(factor.name)),control = list(fnscale=1))
  }
  if (is.character(target.optim)){
    retour.x=as.vector(opti$par) ; names(retour.x)=factor.name ; retour.x=round(retour.x,3)
    retour.y=opti$value
  }else{
    retour.x=as.vector(opti$par) ; names(retour.x)=factor.name ; retour.x=round(retour.x,3)
    vec=as.vector(opti$par)
    pseudo.D=as.data.frame(t(as.matrix(vec))) ; colnames(pseudo.D)=factor.name ;pseudo.formula=as.formula(paste(as.character(formula)[1],as.character(formula)[3],sep=""))
    pseudo.mat=model.matrix(pseudo.formula,pseudo.D)
    y=as.numeric(as.matrix(pseudo.mat)%*%b[colnames(pseudo.mat)])
    retour.y=y
  }
  retour=list(retour.x,retour.y)
  names(retour)=c("factor.values","response.value")
  return(retour)
}
