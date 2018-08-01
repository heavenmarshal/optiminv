evallogpost <- function(y,xi,tlen)
{
    dev <- drop(crossprod(y-xi))/tlen
    loglik <- -0.5*log(2*pi)-0.5*tlen
    loglik <- loglik - 0.5 * tlen* log(dev)
}

lasvdgpEsl2d <- function(design, resp, xi, nstarts=1, n0=5, nn=10,
                         nfea=min(1000, nrow(design)), resvdThres=min(5, nn-n0),
                         every=min(5, nn-n0), frac=.95, gstart=1e-4,
                         maxit=10000, tmax=10, ti=10, kersig=1.0/ti, nthread=1)
{
    ndesign <- nrow(design)
    nparam <- ncol(design)
    tlen <- ncol(resp)
    nsamples <- maxit+1
    logpost <- apply(resp,2,evallogpost,xi,tlen)
    post <- exp(logpost)
    idx <- sample(1:ndesign,nstarts,replace=TRUE,prob=post)
    xstarts <- design[idx,,drop=FALSE]
    out <- .C("lasvdgpEsl2d", as.integer(nparam), as.integer(maxit),
              as.integer(ndesign), as.integer(tlen), as.integer(n0),
              as.integer(nn), as.integer(nfea), as.integer(resvdThres),
              as.integer(every), as.integer(nthread), nas.integer(nstarts),
              as.double(frac), as.double(gstart), as.double(ti),
              as.double(kersig), as.double(t(xstarts)), as.double(xi),
              as.double(t(design)), as.double(resp), mins = double(nstarts),
              argmins=double(nstarts*nparam), fvals = double(nstarts*nsamples),
              samples=double(nstarts*nparam*nsamples))
    fvals <- matrix(out$fvals,nrow=nsamples, ncol=nstarts)
    samples <- array(out$samples, c(nparam,nsamples,nstarts))
    ret <- list(mins=out$mins, argmins=out$argmins, fvals = fvals,
                samples = samples)
    return(ret)
}
lasvdgpEI <- function(design, resp, xi, n0=5, nn=10, nfea=min(1000, nrow(design)),
                      resvdThres=min(5, nn-n0), every=min(5, nn-n0), frac=.95,
                      gstart=1e-4, maxit=10000, tmax=10, ti=10, kersig=1.0/ti)
{
    ndesign <- nrow(design)
    nparam <- ncol(design)
    tlen <- ncol(resp)
    nsamples <- maxit+1
    logpost <- apply(resp,2,evallogpost,xi,tlen)
    post <- exp(logpost)
    idx <- sample(1:ndesign,1,replace=TRUE,prob=post)
    xstart <- design[idx,,drop=FALSE]
    kmin <- min(apply((xi-resp)^2,2,sum))
    out <- .C("lasvdgpEI", as.integer(nparam), as.integer(maxit), as.integer(tmax),
              as.integer(ndesign), as.integer(tlen), as.integer(n0), as.integer(nn),
              as.integer(nfea), as.integer(resvdThres), as.integer(every),
              as.double(frac), as.double(gstart), as.double(ti), as.double(kersig),
              as.double(xstart), as.double(kmin), as.double(xi), as.double(t(design)),
              as.double(resp), mins=double(1), argmin = double(nparam),
              fvals=double(nsamples), samples = double(nparam*nsamples))
    samples <- matrix(out$samples,nrow=nsamples,byrow=TRUE)
    ret <- list(maxei = -out$mins, argmax=out$argmin, fvals=out$fvals,
                samples = samples)
    return(ret)
}
evalLasvdgpEI <- function(design, resp, candid, xi, nn=10, nfea=min(1000, nrow(design)),
                          resvdThres=min(5, nn-n0), every=min(5, nn-n0), frac=.95,
                          gstart=1e-4, nthread=4)
{
    ndesign <- nrow(design)
    ncand <- nrow(candid)
    nparam <- ncol(design)
    tlen <- ncol(resp)
    kmin <- min(apply((xi-resp)^2,2,sum))
    out <- .C("evalLasvdgpEI", as.integer(nparam), as.integer(tlen),
              as.integer(ncand), as.integer(n0), as.integer(nn),
              as.integer(nfea), as.integer(resvdThres), as.integer(every),
              as.integer(nthread), as.double(frac), as.double(gstart),
              as.double(kmin), as.double(xi), as.double(t(design)),
              as.double(resp), as.double(t(candid)), eival=double(ncand))
    return(out$eival)
}
