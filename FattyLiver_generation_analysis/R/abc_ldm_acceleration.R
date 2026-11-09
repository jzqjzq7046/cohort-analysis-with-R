# Algebraically exact acceleration of the published matched-set LDM statistic.
# Within-set-constant nuisance columns have zero crossproduct with the already
# set-adjusted response, and remain constant under the permitted permutations.
# The official statistic is still called on the remaining columns. No test,
# residual, permutation, denominator, transformation or inference is changed.
abc_ldm_acceleration <- function(family,run,verify_every=1000L) {
  ns<-asNamespace('LDM');original<-get('ldm.stat',ns);original_ldm<-get('ldm',ns);family<-as.character(family)
  stopifnot(as.character(packageVersion('LDM'))=='6.0.1')
  original_x<-old_low<-old_up<-keys<-key_column<-NULL;calls<-0L;checks<-list()
  reduce<-function(x,low,up){
      stopifnot(length(low)==2L,is.null(original_x));original_x<<-x;old_low<<-low;old_up<<-up
      nuisance<-x[,seq.int(low[1],up[1]),drop=FALSE]
      stopifnot(max(abs(crossprod(nuisance)-diag(ncol(nuisance))))<1e-8)
      means<-rowsum(nuisance,family,reorder=FALSE);n<-as.numeric(table(factor(family,levels=rownames(means))))
      means<-means/n;within<-nuisance-means[match(family,rownames(means)),,drop=FALSE]
      dec<-svd(within);k<-sum(dec$d>1e-8);stopifnot(k>0,k<ncol(nuisance))
      # Orthogonal rotation isolates the family-invariant null directions. Their
      # norm contribution is zero after the package's nuisance projection.
      fastx<-cbind(dec$u[,seq_len(k),drop=FALSE],x[,seq.int(low[2],up[2]),drop=FALSE])
      key_column<<-which(vapply(seq_len(ncol(fastx)),function(i)!anyDuplicated(fastx[,i]),logical(1)))[1]
      stopifnot(is.finite(key_column));keys<<-fastx[,key_column]
      list(x=fastx,low=c(1L,k+1L),up=c(k,ncol(fastx)))
  }
  accelerated<-function(x,low,up,resid,ss.tot,adjust.for.confounders,comp.anal=FALSE,comp.anal.adjust='median',comp.effect=NULL){
    calls<<-calls+1L;stopifnot(adjust.for.confounders,!comp.anal)
    fast<-original(x=x,low=low,up=up,resid=resid,ss.tot=ss.tot,
      adjust.for.confounders=adjust.for.confounders,comp.anal=comp.anal,comp.anal.adjust=comp.anal.adjust,comp.effect=comp.effect)
    if(calls<=4L||calls%%verify_every==0L){
      permutation<-match(x[,key_column],keys);stopifnot(!anyNA(permutation),!anyDuplicated(permutation))
      ref<-original(x=original_x[permutation,,drop=FALSE],low=old_low,up=old_up,resid=resid,ss.tot=ss.tot,adjust.for.confounders=adjust.for.confounders,
        comp.anal=comp.anal,comp.anal.adjust=comp.anal.adjust,comp.effect=comp.effect)
      err<-max(abs(unlist(fast)-unlist(ref)));scale<-max(1,max(abs(unlist(ref))))
      checks[[length(checks)+1L]]<<-data.frame(statistic_call=calls,original_columns=ncol(original_x),retained_columns=ncol(x),
        maximum_absolute_difference=err,relative_error=err/scale)
      if(!is.finite(err)||err/scale>1e-10)stop('LDM acceleration failed equivalence check')
    }
    fast
  }
  # Insert the algebraic basis reduction once, after official response
  # projection, before observed statistics and the unmodified permutation loop.
  txt<-deparse(body(original_ldm),width.cutoff=500);at<-grep('ldm.obs.freq = ldm.stat(',txt,fixed=TRUE);stopifnot(length(at)==1L)
  insertion<-c('abc_reduced <- abc_reduce_design(x.design, low, up)','x.design <- abc_reduced$x','low <- abc_reduced$low','up <- abc_reduced$up')
  local_ldm<-original_ldm;body(local_ldm)<-parse(text=append(txt,insertion,after=at-1L))[[1]]
  env<-new.env(parent=environment(original_ldm));env$abc_reduce_design<-reduce;environment(local_ldm)<-env
  replace<-function(name,value){unlockBinding(name,ns);assign(name,value,ns);lockBinding(name,ns)}
  replace('ldm.stat',accelerated);replace('ldm',local_ldm)
  on.exit({replace('ldm.stat',original);replace('ldm',original_ldm)},add=TRUE)
  result<-run()
  list(result=result,equivalence=do.call(rbind,checks),calls=calls,official_statistic_hash=digest::digest(body(original),algo='sha256'),
    official_ldm_hash=digest::digest(body(original_ldm),algo='sha256'))
}
