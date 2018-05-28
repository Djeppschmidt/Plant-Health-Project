\name{diffDESeq}
\alias{diffDESeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{diffDESeq help
%%  ~~function to do ... ~~
}
\description{
%%  ~~ conducts DESeq2 variance stabilization for PHP analysis. ~~
}
\usage{
diffDESeq(x)
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
diffDESeq<-function(ps){
  L1 <-phyloseq_to_deseq2(ps, ~ group + Sampling_date + group:Sampling_date)

  diagdds = DESeq(L1, test="Wald", fitType="parametric")

  res = results(diagdds)

  return(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
