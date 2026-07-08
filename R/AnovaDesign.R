#' ANOVA tables of \pkg{DoETeachR} objects
#'
#' Computes the (type III) ANOVA table of a \pkg{DoETeachR} object using the model from \bold{formula} and the corresponding adjusted R-squared.
#'
#' @param formula A formula that defines the model. Only order two effects are supported at maximum.
#' @param design A design from \code{\link[DoETeachR]{gen.FrFD}}, \code{\link[DoETeachR]{gen.FFD}} or \code{\link[DoETeachR]{gen.CCD}}.
#' @param alpha.risk.model The alpha risk of the global test of the model. If the global test shows a p.value larger than \bold{alpha.risk.model}, then all effects from the model are considered non significant.
#'
#' @returns A list where the first element is the ANOVA table and the second is the adjusted R-squared of the model.
#'
#' @export
#'
#' @examples
#'
#' design=gen.FFD(responses=c("R1"),factors=LETTERS[1:4])
#' design$R1=runif(nrow(design))
#' AnovaDesign(R1~.*.,design)

AnovaDesign=function(formula,design,alpha.risk.model=0.05){
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
  if (is.numeric(alpha.risk.model) | is.integer(alpha.risk.model)){
    if (length(alpha.risk.model)==1){
      if (alpha.risk.model>1 | alpha.risk.model<0){
        stop("alpha.risk.model must be between 0 and 1")
      }
    }else{
      stop("length(alpha.risk.model) must equal 1")
    }
  }else{
    stop("class(alpha.risk.model) must be integer or numeric")
  }
  mod=lm(formula,design) ; x=model.matrix(mod) ; b=coef(mod)
  if (any(sapply(strsplit(colnames(x),"[:]"),length)>2)){
    stop("formula must no contain effect of order larger than 2")
  }
  if (mod$df.residual==0){
    stop("No degree of freedom for the error estimation")
  }
  factor.name=colnames(design)[attr(design,"colonne")=="factor"] ; factor.name=factor.name[apply(sapply(factor.name, gregexpr,text=colnames(x))>0,2,any)]
  factor.pure=interaction(design[,factor.name]) ; y.pure=design[,as.character(formula)[2]]
  mod.pure=lm(y.pure~factor.pure)
  ali=alias(mod)
  if (!is.null(ali$Complete)){
    ali=round(as.matrix(as.data.frame(ali$Complete)),12)
    for (j in 1:ncol(ali)){
      if (any(ali[,j]!=0)){
        ou.1=which(ali[,j]==1)
        premier=colnames(ali)[j]
        seconds=rownames(ali)[ou.1]
        names(b)[names(b)==premier]=paste(premier,paste(seconds,collapse = "+"),sep = "+")
        b=b[-which(names(b)%in%seconds)]
        colnames(x)[colnames(x)==premier]=paste(premier,paste(seconds,collapse = "+"),sep = "+")
        x=x[,-which(colnames(x)%in%seconds),drop=FALSE]
      }
    }
  }
  if (mod.pure$df.residual!=0 & ncol(x)>2){
    SS.erreur.pure=as.numeric(crossprod(mod.pure$residuals)) ; df.erreur.pure=mod.pure$df.residual
    fit=x[,-1,drop=FALSE]%*%as.matrix(b[-1]) ; df.fit=ncol(x[,-1,drop=FALSE])
    SS.fit=as.numeric(crossprod(fit-x[,1,drop=FALSE]%*%solve(crossprod(x[,1,drop=FALSE]))%*%crossprod(x[,1,drop=FALSE],fit)))
    p.fit=pf((SS.fit/df.fit)/(SS.erreur.pure/df.erreur.pure),df.fit,df.erreur.pure,lower.tail = FALSE)
    if (p.fit<=alpha.risk.model){
      retour.an=matrix(NA,ncol(x)+1,5) ; rownames(retour.an)=c(colnames(x)[-1],"Lack-Of-Fit","Pure Error") ; colnames(retour.an)=c("SSQ","Df","MS","F","p.value")
      for (i in 1:(nrow(retour.an)-2)){
        nom=rownames(retour.an)[i]
        xo=x[,!colnames(x)%in%nom,drop=FALSE]
        fiti=x[,nom,drop=FALSE]*b[nom]
        SSi=crossprod(fiti-xo%*%solve(crossprod(xo))%*%crossprod(xo,fiti)) ; dfi=ncol(x[,nom,drop=FALSE])
        retour.an[nom,1]=SSi ; retour.an[nom,2]=dfi
      }
      SS.erreur=as.numeric(crossprod(mod$residuals)) ; df.erreur=mod$df.residual
      retour.an[nrow(retour.an)-1,1]=SS.erreur-SS.erreur.pure ; retour.an[nrow(retour.an)-1,2]=df.erreur-df.erreur.pure
      retour.an[nrow(retour.an),1]=SS.erreur.pure ; retour.an[nrow(retour.an),2]=df.erreur.pure
      retour.an[,3]=retour.an[,1]/retour.an[,2] ; retour.an[,4]=retour.an[,3]/retour.an[nrow(retour.an),3]
      retour.an[,5]=pf(retour.an[,4],retour.an[,2],retour.an[nrow(retour.an),2],lower.tail = FALSE)
      retour.an[,1:4]=round(retour.an[,1:4],5)
      retour.an[nrow(retour.an),4:5]=""
      retour.an=as.data.frame(retour.an)
      etoiles = rep("   ", length(retour.an$p.value))
      retour.an$p.value=as.numeric(retour.an$p.value)
      etoiles[retour.an$p.value <= 0.10] = ".  "
      etoiles[retour.an$p.value <= 0.05] = "*  "
      etoiles[retour.an$p.value <= 0.01] = "** "
      etoiles[retour.an$p.value <= 0.001] = "***"
      retour.an$Sig=etoiles
      names(retour.an)[6]="" ; retour.an[nrow(retour.an),6]=""
      retour.an$p.value=format.pval(as.numeric(retour.an$p.value),na.form = "")
      SST=as.numeric(crossprod(y.pure-mean(y.pure)))
      if (df.erreur>0){
        retour.R2.adj=max(0,1-((SS.erreur/df.erreur)/(SST/(length(y.pure)-1))))
        retour.R2.adj=round(retour.R2.adj,3)
      }else{
        retour.R2.adj=0
      }
    }else{
      retour.an=retour.R2.adj=paste("Non significant model, p.value model = ",format.pval(p.fit),sep="")
    }
    return(list(Anova=retour.an,R2.adj=retour.R2.adj))
  }else{
    SS.erreur=as.numeric(crossprod(mod$residuals)) ; df.erreur=mod$df.residual
    fit=x[,-1,drop=FALSE]%*%as.matrix(b[-1]) ; df.fit=ncol(x[,-1,drop=FALSE])
    SS.fit=as.numeric(crossprod(fit-x[,1,drop=FALSE]%*%solve(crossprod(x[,1,drop=FALSE]))%*%crossprod(x[,1,drop=FALSE],fit)))
    p.fit=pf((SS.fit/df.fit)/(SS.erreur/df.erreur),df.fit,df.erreur,lower.tail = FALSE)
    if (p.fit<=alpha.risk.model){
      retour.an=matrix(NA,ncol(x),5) ; rownames(retour.an)=c(colnames(x)[-1],"Total Error") ; colnames(retour.an)=c("SSQ","Df","MS","F","p.value")
      for (i in 1:(nrow(retour.an)-1)){
        nom=rownames(retour.an)[i]
        xo=x[,!colnames(x)%in%nom,drop=FALSE]
        fiti=x[,nom,drop=FALSE]*b[nom]
        SSi=crossprod(fiti-xo%*%solve(crossprod(xo))%*%crossprod(xo,fiti)) ; dfi=ncol(x[,nom,drop=FALSE])
        retour.an[nom,1]=SSi ; retour.an[nom,2]=dfi
      }
      retour.an[nrow(retour.an),1]=SS.erreur ; retour.an[nrow(retour.an),2]=df.erreur
      retour.an[,3]=retour.an[,1]/retour.an[,2] ; retour.an[,4]=retour.an[,3]/retour.an[nrow(retour.an),3]
      retour.an[,5]=pf(retour.an[,4],retour.an[,2],retour.an[nrow(retour.an),2],lower.tail = FALSE)
      retour.an[,1:4]=round(retour.an[,1:4],5)
      retour.an[nrow(retour.an),4:5]=""
      retour.an=as.data.frame(retour.an)
      etoiles = rep("   ", length(retour.an$p.value))
      retour.an$p.value=as.numeric(retour.an$p.value)
      etoiles[retour.an$p.value <= 0.10] = ".  "
      etoiles[retour.an$p.value <= 0.05] = "*  "
      etoiles[retour.an$p.value <= 0.01] = "** "
      etoiles[retour.an$p.value <= 0.001] = "***"
      retour.an$Sig=etoiles
      names(retour.an)[6]="" ; retour.an[nrow(retour.an),6]=""
      retour.an$p.value=format.pval(as.numeric(retour.an$p.value),na.form = "")
      SST=as.numeric(crossprod(y.pure-mean(y.pure)))
      if (df.erreur>0){
        retour.R2.adj=max(0,1-((SS.erreur/df.erreur)/(SST/(length(y.pure)-1))))
        retour.R2.adj=round(retour.R2.adj,3)
      }else{
        retour.R2.adj=0
      }
    }else{
      retour.an=retour.R2.adj=paste("Non significant model, p.value model = ",format.pval(p.fit),sep="")
    }
    return(list(Anova=retour.an,R2.adj=retour.R2.adj))
  }
}
