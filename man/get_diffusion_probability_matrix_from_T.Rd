% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_diff_prob_matrix.R
\name{get_diffusion_probability_matrix_from_T}
\alias{get_diffusion_probability_matrix_from_T}
\alias{get_diffu_Pt_from_T}
\alias{get_diffu_Pt_from_Pi}
\alias{get_diffusion_probability_matrix_from_Pi}
\title{Diffusion probability matrix from transition matrix}
\usage{
get_diffusion_probability_matrix_from_T(Pi, tau)

get_diffu_Pt_from_T(Pi, tau)

get_diffu_Pt_from_Pi(Pi, tau)

get_diffusion_probability_matrix_from_Pi(Pi, tau)
}
\arguments{
\item{Pi}{Transition matrix. We do not use \code{T} to avoid conflicts
with the abbreviation for \code{TRUE}), instead we indicated the transition
matrix with the capital greek letter \eqn{\Pi} in the equations and
\code{Pi} in the code.}

\item{tau}{diffusion time}
}
\value{
\eqn{exp^{-\tau (I - \Pi)}}, exponential of the normalized Laplacian
matrix corresponding to the given transition rate matrix (or transition
probability matrix of a discrete-time Markov chain).
}
\description{
Description here
}
\references{
De Domenico, M. (2017). Diffusion Geometry Unravels the Emergence of
Functional Clusters in Collective Phenomena. Physical Review Letters.
\doi{10.1103/PhysRevLett.118.168301}

Bertagnolli, G., & De Domenico, M. (2020). Diffusion Geometry of Multiplex
and Interdependent Systems.
\href{https://arxiv.org/abs/2006.13032}{arxiv preprint arxiv:2006.13032}
}
\seealso{
\code{\link{get_diffusion_probability_matrix}}
}
\keyword{Markov-chain}
\keyword{probabilities}
\keyword{transition}
