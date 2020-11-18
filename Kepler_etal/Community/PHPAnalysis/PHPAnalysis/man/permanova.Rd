\name{permanova}
\alias{permanova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{permanova help
%%  ~~function to do ... ~~
}
\description{
%%  ~~ runs permaova for PHP analysis, and subsets the results objects for downstream analysis. ~~
}
\usage{
permanova(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
permanova<-function(ps, n, indics) {
  require(phyloseq)
  require(vegan)
  require(stats)
  out=adonis(otu_table(ps) ~ System.loc * Glyphosphate_Treatment * Sampling_date * Soil_Zone, strata = sample_data(ps)$Loc_plot_ID, as(sample_data(ps), "data.frame"))

  tcoeffs<-data.frame(t(out$coefficients))
  tcoeffs$ID<-rownames(tcoeffs)

  syst=rbind(tcoeffs[which(tcoeffs$System.loc1> quantile(tcoeffs$System.loc1, prob=1-n/100)),], tcoeffs[which(tcoeffs$System.loc1< quantile(tcoeffs$System.loc1, prob=n/100)),]) #filter taxa based on steepest n% slopes in system.loc

  glyph=rbind(tcoeffs[which(tcoeffs$Glyphosphate_Treatment1> quantile(tcoeffs$Glyphosphate_Treatment1, prob=1-n/100)),], tcoeffs[which(tcoeffs$Glyphosphate_Treatment1< quantile(tcoeffs$Glyphosphate_Treatment1, prob=n/100)),]) #filter taxa based on steepest 5% slopes in Glyphosphate Treatment

  date=rbind(tcoeffs[which(tcoeffs$Sampling_date1> quantile(tcoeffs$Sampling_date1, prob=1-n/100)),], tcoeffs[which(tcoeffs$Sampling_date1< quantile(tcoeffs$Sampling_date1, prob=n/100)),]) #filter taxa based on steepest 5% slopes in Sampling Date

  a<-NULL
  b<-NULL
  c<-NULL
  d<-NULL
  e<-NULL
  ind<-as.data.frame(a,b,c,d,e)

  for(i in indics){
    first<-tcoeffs[which(tcoeffs$ID==i),]
    ind<-rbind(ind, first)
  }
  rownames(ind)<-indics # testing...
  ind<-data.frame(ind, factor(indics)) #not sure if this will work
  return(list("out"=out, "syst"=syst, "glyph"=glyph, "date"=date, "sppInt"=ind))

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
