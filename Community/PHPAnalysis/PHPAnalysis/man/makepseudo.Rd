\name{makepseudo}
\alias{makepseudo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{makepseudo help
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makePseudo(ps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ps}{
%%
}
}
\details{
%%  ~~ phyloseq object should be of bacterial organisms, and contain taxonomy assignments down to the genus level ~~
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
  udo<-subset_taxa(ps, Genus=="Pseudomonas")
  udo2<-as.matrix(tax_table(udo))
  indics<-rownames(udo2)
  indics
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
