\name{makefusarium}
\alias{makefusarium}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{makefusarium help
%%  ~~function to do ... ~~
}
\description{
%%  ~~ takes a phyloseq object containing a bacterial community, and returns a subset containing all pseudomonas as a phyloseq object ~~
}
\usage{
makefusarium(x)
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
{
  require(phyloseq)
  udo<-subset_taxa(ps, Genus=="g__Fusarium")
  udo2<-as.matrix(tax_table(udo))
  indics<-rownames(udo2)
  indics
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
