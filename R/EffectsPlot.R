#' Plot Effects of \pkg{DoETeachR} objects
#'
#' Plots either direct or two-way interactions effects from \bold{formula} of a \pkg{DoETeachR} object
#'
#' @param formula A formula that defines the model. Only order two effects are supported at maximum.
#' @param design A design from \code{\link[DoETeachR]{gen.FrFD}}, \code{\link[DoETeachR]{gen.FFD}} or \code{\link[DoETeachR]{gen.CCD}}.
#' @param effects Either `"direct"` or `"interaction"` specifying if direct or interaction effects must be plotted.
#' @param bounds Either NULL (Default) or a numeric vector of length two defining the range of values that must be considered for the factors. If NULL, the range of values for the factors are the minimum and maximum values of the plotted factors present in \bold{formula}.
#'
#' @returns An effect plot
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'design=gen.FFD(responses=c("R1"),factors=LETTERS[1:4])
#'design$R1=runif(nrow(design))
#'EffectsPlot(R1~.*.,design,"direct")
#'EffectsPlot(R1~.*.,design,"interaction")


EffectsPlot=function(formula,design,effects="direct",bounds=NULL){
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
  if (is.character(effects)){
    if (length(effects)==1){
      if (!effects%in%c("direct","interaction")){
        stop("effects must be direct or interaction")
      }
    }else{
      stop("length(effects) must be one")
    }
  }else{
    stop("class(effects) must be character")
  }
  mod=lm(formula,design) ; x=model.matrix(mod) ; b=coef(mod)
  verif.alias=alias(mod)
  if (!is.null(verif.alias$Complete)){
    stop("There are aliased term in the model")
  }
  if (any(sapply(strsplit(colnames(x),"[:]"),length)>2)){
    stop("formula must no contain effect of order larger than 2")
  }
  if (is.numeric(bounds) | is.null(bounds)){
    if (!is.null(bounds)){
      if (length(bounds)!=2){
        stop("length(bounds) must be 2 when class(bounds) is numeric")
      }
    }
  }else{
    stop("class(bounds) must be numeric or NULL")
  }
  if (is.null(bounds)){
    bounds=c(min(x),-min(x))
  }
  factor.name=colnames(design)[attr(design,"colonne")=="factor"] ; factor.name=factor.name[apply(sapply(factor.name, gregexpr,text=colnames(x))>0,2,any)]
  pseudo.formula=as.formula(paste(as.character(formula)[1],as.character(formula)[3],sep=""))
  ou.inter=sapply((gregexpr("[:]",colnames(x))),function(x) x>0)
  if (effects=="direct"){
    grille=matrix(rep(0,500*length(factor.name)),nrow = 500) ; colnames(grille)=factor.name
    df.plot=NULL
    for (f in factor.name){
      grille.f=grille ; grille.f[,f]=seq(bounds[1],bounds[2],length.out=500) ; x.f=model.matrix(pseudo.formula,as.data.frame(grille.f))
      df.plot=rbind(df.plot,cbind(f,x.f[,f],x.f%*%b[colnames(x.f)]))
    }
    colnames(df.plot)=c("Factor","Value","Response") ; df.plot=as.data.frame(df.plot) ; df.plot[,1]=factor(df.plot[,1],levels = unique(df.plot[,1])) ; df.plot[,2]=as.numeric(df.plot[,2]) ; df.plot[,3]=as.numeric(df.plot[,3])
    g=ggplot(df.plot,mapping = aes(x=df.plot[,2],y=df.plot[,3]))+theme_bw()
    g=g+geom_line(linewidth=1.2,col="#9400D3")+facet_wrap(~Factor,ncol=nlevels(df.plot[,1]),strip.position="bottom") + ylab(colnames(mod$model)[1])
    g = g + theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 12, face = "bold"),axis.text = element_text(size=max(5,14-length(factor.name))),strip.text = element_text(size = 12, face = "bold"),strip.background = element_rect(fill = "deepskyblue1"))
    g=g+geom_hline(yintercept = b[1],linetype="dotted",linewidth=1)
  }else{
    interactions=colnames(x[,ou.inter,drop=FALSE])
    grille=matrix(rep(0,1000*length(factor.name)),nrow = 1000) ; colnames(grille)=factor.name
    df.plot=NULL
    for (i in interactions){
      dans.i=strsplit(i,"[:]")[[1]]
      grille.i=grille ; grille.i[,dans.i[1]]=rep(seq(bounds[1],bounds[2],length.out=500),2) ; grille.i[,dans.i[2]]=rep(c(bounds[1],bounds[2]),each=500)
      x.i=model.matrix(pseudo.formula,as.data.frame(grille.i))
      df.plot=rbind(df.plot,cbind(i,x.i[,dans.i[1]],x.i[,dans.i[2]],x.i%*%b[colnames(x.i)]))
    }
    colnames(df.plot)=c("Interaction","Value1","Value2","Response") ; df.plot=as.data.frame(df.plot) ; df.plot[,1]=factor(df.plot[,1],levels = unique(df.plot[,1])) ; df.plot[,2]=as.numeric(df.plot[,2]) ; df.plot[,3]=as.factor(round(as.numeric(df.plot[,3]),3)) ; df.plot[,4]=as.numeric(df.plot[,4])
    g=ggplot(df.plot,mapping = aes(x=df.plot[,2],y=df.plot[,4],group = df.plot[,3],color=df.plot[,3]))+theme_bw()
    g=g+geom_line(linewidth=1.2)+facet_wrap(~Interaction,ncol=nlevels(df.plot[,1]),strip.position="bottom") + ylab(colnames(mod$model)[1])
    g = g + theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 12, face = "bold"),axis.text = element_text(size=max(5,14-length(interactions))),strip.text = element_text(size = 12, face = "bold"),strip.background = element_rect(fill = "deepskyblue1"),legend.position = "top",legend.title = element_text(size = 12, face = "bold"),legend.text = element_text(face="bold",size = 12))
    g = g + scale_color_manual(name="Right Element of the Interaction:",values=c("blue","red"))
    g=g+geom_hline(yintercept = b[1],linetype="dotted",linewidth=1)
  }
  return(g)
}
