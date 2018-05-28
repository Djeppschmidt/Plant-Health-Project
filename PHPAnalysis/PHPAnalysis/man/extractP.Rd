\name{extractP}
\alias{extractP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{extractP help
%%  ~~function to do ... ~~
}
\description{
%%  ~~ extracts p values from DESeq2 logfold change results object ~~
}
\usage{
extractP(x)
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
extractP<-function(ps, pvalue){
  if (identical(rownames(tax_table(ps)), rownames(pvalue))) {
    T1<-ps
    new_tax<-as.data.frame(as.matrix(tax_table(T1)))
    new_tax$DESeq_padj<-pvalue$padj
    tax_table(T1)<-tax_table(as.matrix(new_tax))
    T1
  }
  else {return("rownames out of order")}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
