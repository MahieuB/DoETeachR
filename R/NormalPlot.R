#' Plot Normality of Effects of \pkg{DoETeachR} objects
#'
#' Plots the empirical cumulative probabilities against the absolute effects from \bold{formula} of a \pkg{DoETeachR} object. Adds a reference line following Zahn, D. (1975) method.
#'
#'
#' @param formula A formula that defines the model. Only order two effects are supported at maximum.
#' @param design A design from \code{\link[DoETeachR]{gen.FrFD}} that has no degree of freedom for the residuals given the model of \bold{formula}
#'
#' @returns A plot of normality of effects.
#'
#' @import ggrepel
#' @import ggplot2
#' @import stats
#' @import MASS
#'
#' @references Zahn, D (1975) Modifications of and Revised Critical Values for the Half-Normal Plot. Technometrics 17(2), 189-200
#'
#' @export
#'
#' @examples
#' design=gen.FrFD(responses=c("R1"),factors=LETTERS[1:4],frac=2^1)
#' design$R1=runif(nrow(design))
#' NormalPlot(R1~.*.,design)

NormalPlot=function(formula,design){
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
  mod=lm(formula,design)
  if (mod$df.residual>0 & is.null(alias(mod)$Complete)){
    stop("df.residual>0 & no coefficent are aliased --> use AnovaDesign")
  }
  b=coef(mod) ; eff=b[-1] ; x=model.matrix(mod)
  if (any(sapply(strsplit(colnames(x),"[:]"),length)>2)){
    stop("formula must no contain effect of order larger than 2")
  }
  ali=stats::alias(mod)
  if (!is.null(ali$Complete)){
    ali=round(as.matrix(as.data.frame(ali$Complete)),12)
    for (j in 1:ncol(ali)){
      if (any(ali[,j]!=0)){
        ou.1=which(ali[,j]==1)
        premier=colnames(ali)[j]
        seconds=rownames(ali)[ou.1]
        names(eff)[names(eff)==premier]=paste(premier,paste(seconds,collapse = "+"),sep = "+")
        eff=eff[-which(names(eff)%in%seconds)]
      }
    }
  }
  empi=(rank(eff)-0.5)/length(eff)
  df.plot=data.frame(eff,empi)
  g=ggplot(df.plot,mapping=aes(x=df.plot[,1],y=df.plot[,2]))+theme_bw()+ylim(0,1)
  combien.droite=floor(0.683 * length(eff))
  choix.droite=rank(abs(eff))<=combien.droite
  mod.droite=lm(empi[choix.droite]~eff[choix.droite])
  g=g+geom_point(size=3)+geom_abline(intercept = coef(mod.droite)[1],slope=coef(mod.droite)[2],linewidth=1,col="blue")
  g = g + theme(axis.title = element_text(size = 14, face = "bold"),axis.text = element_text(size=12),strip.text = element_text(size = 12, face = "bold"),strip.background = element_rect(fill = "deepskyblue1"),legend.position = "top",legend.title = element_text(size = 12, face = "bold"),legend.text = element_text(face="bold",size = 12))
  g = g +xlab("Effect") + ylab("Empirical Cumulative Probabilities")
  lab=df.plot
  g=g+geom_label_repel(as.data.frame(lab),mapping=aes(x=lab[,1],y=lab[,2],label=rownames(lab)),linewidth=0,colour="black",size=5,segment.size=1,label.padding = 0,
                       min.segment.length = 1,nudge_y = -0.025)
  return(g)
}
