# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath)
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
loadNamespace("ANCOMBC")
stopifnot(as.character(packageVersion("ANCOMBC")) == "2.12.0")
replacement <- new.env(parent=asNamespace("ANCOMBC"))
sys.source(args[[2]], envir=replacement)
assignInNamespace(".bias_em",replacement$.bias_em, ns="ANCOMBC")
results <- list()
for (replicate in seq_len(32)) {
  reads <- as.matrix(read.delim(file.path(input,paste0("replicate-",replicate,".tsv")),row.names=1,check.names=FALSE))
  meta <- data.frame(condition=factor(c(rep("reference",4),rep("contrast",16)),levels=c("reference","contrast")),
                     unused=seq_len(20),row.names=colnames(reads))
  diagnostics <- list()
  trace(".bias_em",where=asNamespace("ANCOMBC"),print=FALSE,exit=quote({
    .GlobalEnv$diagnostics[[length(.GlobalEnv$diagnostics)+1L]] <- c(iterations=iterNum,residual=epsilon)
  }))
  fit <- tryCatch(ANCOMBC::ancombc(data=reads,aggregate_data=reads,meta_data=meta,
    tax_level=NULL,formula="condition",p_adj_method="BH",pseudo=1,prv_cut=.1,lib_cut=0,
    group="condition",struc_zero=FALSE,neg_lb=FALSE,tol=1e-7,max_iter=20000,
    conserve=TRUE,global=FALSE,n_cl=1,verbose=FALSE),error=function(e)e)
  untrace(".bias_em",where=asNamespace("ANCOMBC"))
  if(inherits(fit,"error")) {
    results[[replicate]] <- data.frame(replicate,status=conditionMessage(fit),false_calls=NA,iterations=NA,residual=NA,native_false_calls=NA,max_q_gap=NA)
  } else {
    diag <- do.call(rbind,diagnostics)
    converged <- nrow(diag)>0 && all(diag[,"iterations"]<20000) && all(diag[,"residual"]<=1e-7)
    q <- fit$res$q_val$conditioncontrast
    native_file <- file.path(input,paste0("native-",replicate,".csv"))
    native_calls <- max_q_gap <- NA_real_
    if (file.exists(native_file)) {
      native <- read.csv(native_file,colClasses=c(feature="character"))
      ids <- if("taxon" %in% names(fit$res$q_val)) as.character(fit$res$q_val$taxon) else rownames(fit$res$q_val)
      stopifnot(setequal(ids,native$feature),all(is.finite(q)))
      native_calls <- sum(native$qvalue<=.05)
      max_q_gap <- max(abs(native$qvalue-q[match(native$feature,ids)]))
    }
    results[[replicate]] <- data.frame(replicate,status=if(converged) "converged" else "unconverged",
       false_calls=sum(q<=.05,na.rm=TRUE),iterations=max(diag[,"iterations"]),residual=max(diag[,"residual"]),native_false_calls=native_calls,max_q_gap)
  }
  write.csv(do.call(rbind,results),file.path(dirname(input),"small-unbalanced-repaired-R.csv"),row.names=FALSE)
  cat("Finished replicate",replicate,"\n"); flush.console()
}
