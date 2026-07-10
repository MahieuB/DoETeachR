#' Plot Response Surface of \pkg{DoETeachR} objects
#'
#' Plots the response surface from \bold{formula} as a function of \bold{plotted.factors} of a \pkg{DoETeachR} object. A green point is added at the optimal value depending on the target defined by \bold{target.optim}.
#'
#' @param formula A formula that defines the model. Only order two effects are supported at maximum.
#' @param design A design from \code{\link[DoETeachR]{gen.FrFD}}, \code{\link[DoETeachR]{gen.FFD}} or \code{\link[DoETeachR]{gen.CCD}}.
#' @param plotted.factors Either NULL (Default) of character vector of length two naming the factors of \bold{formula} for which the response surface must be plotted. If NULL, the two plotted factors are those whose sum of squares of effects are the largest.
#' @param bounds Either NULL (Default) or a numeric vector of length two defining the range of values that must be considered for the \bold{plotted.factors}. If NULL, the range of values are the minimum and maximum values of the plotted factors present in \bold{formula}.
#' @param target.optim Either `"min"`, `"max"` or a numeric value corresponding to a target value of the response indicating whether the objective is to minimize, maximize or maintain the value of the response at the target respectively.
#' @param value.non.plotted Either `"optimal"` or `"central"` indicating whether the non plotted factors should be considered at their optimal value or their central value when drawing the response surface of \bold{plotted.factors}.
#'
#' @returns A surface response plot
#'
#' @import ggplot2
#' @import ggrepel
#' @import stats
#' @import grDevices
#'
#' @export
#'
#' @examples
#'
#' design=gen.CCD(responses=c("R1"),factors=LETTERS[1:4])
#' design$R1=runif(nrow(design))
#' SurfacePlot(R1~.*.+I(A^2)+I(B^2)+I(C^2)+I(D^2),design)

SurfacePlot=function(formula,design,plotted.factors=NULL,bounds=NULL,target.optim="min",value.non.plotted="optimal"){
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
  if (is.character(plotted.factors) | is.null(plotted.factors)){
    if (!is.null(plotted.factors)){
      if (length(plotted.factors)!=2){
        stop("length(plotted.factors) must be 2 when class(plotted.factors) is character")
      }
      if (any(!plotted.factors%in%factor.name)){
        stop("plotted.factors are not in formula")
      }
    }
  }else{
    stop("class(plotted.factors) must be character or NULL")
  }
  if (is.numeric(bounds) | is.null(bounds)){
    if (!is.null(bounds)){
      if (length(bounds)!=2){
        stop("length(bounds) must be 2 when class(bounds) is numeric")
      }
      if (all(bounds>0)){
        stop("min(bounds) must be negative")
      }
    }
  }else{
    stop("class(bounds) must be numeric or NULL")
  }
  if (is.character(value.non.plotted)){
    if (length(value.non.plotted)==1){
      if (!value.non.plotted%in%c("optimal","central")){
        stop("value.non.plotted must be optimal or central")
      }
    }else{
      stop("length(value.non.plotted) non plotted must be 1")
    }
  }else{
    stop("class(value.non.plotted) must be character")
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
  if (is.null(plotted.factors)){
    SSQ=drop1(mod,~.)[-1,2,drop=FALSE] ; effect.names=rownames(SSQ)
    for (fa in factor.name){
      inornot=sapply(factor.name, gregexpr,text=effect.names)>0
      sm.SSQ=sapply(factor.name, function(x) sum(SSQ[inornot[,x],]))
      plotted.factors=names(sm.SSQ[rank(-sm.SSQ)<=2])
    }
  }
  grille=expand.grid(replicate(2,seq(bounds[1],bounds[2],length.out=500),simplify = FALSE)) ; colnames(grille)=plotted.factors
  if (value.non.plotted=="optimal"){
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
      opti.value=as.vector(opti$par) ; names(opti.value)=factor.name  ; opti.value=round(opti.value,3)
      grille=cbind(grille,as.matrix(rep(1,nrow(grille)))%*%t(as.matrix(opti.value[!names(opti.value)%in%plotted.factors])))
    }else{
      fopt=function(vec){
        pseudo.D=as.data.frame(t(as.matrix(vec))) ; colnames(pseudo.D)=factor.name ;pseudo.formula=as.formula(paste(as.character(formula)[1],as.character(formula)[3],sep=""))
        pseudo.mat=model.matrix(pseudo.formula,pseudo.D)
        y=as.numeric(as.matrix(pseudo.mat)%*%b[colnames(pseudo.mat)])
        return(sqrt((y-target.optim)^2))
      }
      opti=optim(rep(0,length(factor.name)),fopt,method = "L-BFGS-B",lower=rep(bounds[1],length(factor.name)),upper=rep(bounds[2],length(factor.name)),control = list(fnscale=1))
      opti.value=as.vector(opti$par) ; names(opti.value)=factor.name  ; opti.value=round(opti.value,3)
      grille=cbind(grille,as.matrix(rep(1,nrow(grille)))%*%t(as.matrix(opti.value[!names(opti.value)%in%plotted.factors])))
    }
  }else{
    opti.value=rep(0,length(factor.name)) ; names(opti.value)=factor.name  ; opti.value=round(opti.value,3)
    grille=cbind(grille,as.matrix(rep(1,nrow(grille)))%*%t(as.matrix(opti.value[!names(opti.value)%in%plotted.factors])))
  }
  grille=grille[,factor.name]
  x.formula=as.formula(paste(as.character(formula)[1],as.character(formula)[3],sep=""))
  x.grille=model.matrix(x.formula,grille)
  z.grille=x.grille%*%b[colnames(x.grille)]
  df.plot=data.frame(grille[,plotted.factors],z=z.grille)
  grille.rupture=expand.grid(replicate(length(factor.name),seq(bounds[1],bounds[2],length.out=floor(exp(log(500^2)/length(factor.name)))),simplify = FALSE)) ; colnames(grille.rupture)=factor.name
  x.rupture=model.matrix(x.formula,grille.rupture)
  z.rupture=x.rupture%*%b[colnames(x.rupture)]
  rupture=quantile(z.rupture,probs=seq(0,1,by=0.05),type=5)
  color.rupture=colorRampPalette(c("blue","white","red"))(length(rupture)-1)
  color.plot=color.rupture[table(cut(z.grille,breaks = rupture))>0]
  g=ggplot(df.plot,mapping = aes(x=df.plot[,1],y=df.plot[,2],z=df.plot[,3]))+theme_bw()
  g=g+geom_contour_filled(breaks = rupture)+geom_contour(color="black",breaks = rupture)+geom_vline(xintercept = 0,linetype="dashed")+geom_hline(yintercept = 0,linetype="dashed")
  g=g+scale_fill_manual(values=color.plot,name=as.character(formula)[2])
  g=g+scale_x_continuous(expand = c(0,0),limits = bounds*1.05)+scale_y_continuous(expand = c(0,0),limits = bounds*1.05)+coord_fixed()
  g = g + theme(axis.title = element_text(size = 14, face = "bold"),axis.text = element_text(size=12),legend.title = element_text(face="bold",size=12),plot.title = element_text(hjust=0.5,face = "bold",size=14))
  g=g+xlab(plotted.factors[1])+ylab(plotted.factors[2])+ggtitle(paste(paste(names(opti.value)[!names(opti.value)%in%plotted.factors],opti.value[!names(opti.value)%in%plotted.factors],sep = ": "),collapse = " / "))
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
    opti.value=as.vector(opti$par) ; names(opti.value)=factor.name  ; opti.value=round(opti.value,3)
    grille=cbind(grille,as.matrix(rep(1,nrow(grille)))%*%t(as.matrix(opti.value[!names(opti.value)%in%plotted.factors])))
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
    point.optim=as.vector(opti$par) ; names(point.optim)=factor.name ; point.optim=point.optim[plotted.factors]
  }else{
    point.optim=as.vector(opti$par) ; names(point.optim)=factor.name ; point.optim=point.optim[plotted.factors]
  }
  point.optim=as.data.frame(t(as.matrix(point.optim)))
  suppressWarnings({g=g+geom_point(point.optim,mapping=aes(x=point.optim[,1],y=point.optim[,2],z=0),size=3.5,color="green4")})
  suppressWarnings({g=g+geom_label_repel(point.optim,mapping=aes(x=point.optim[,1],y=point.optim[,2],label="Target",z=0),label.size=NA,colour="green4",size=5,segment.size=1,label.padding = 0.15,fill="black",
                                         min.segment.length = 0,nudge_x = -sign(point.optim[,1])*0.125,nudge_y = -sign(point.optim[,2])*0.125)})
  return(g)
}

