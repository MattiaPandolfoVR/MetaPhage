require(Biostrings)
require(parallel)
require(VirFinder)

#example
#inFaFile <- system.file("data", "tara_host.fa", package="VirFinder")
#vir.pred <- parVF.pred(inFaFile, cores=8)


parVF.run <- function (seqFaIn)
{
  seqFa <- strsplit(x=as.character(seqFaIn),split="",fixed=T)[[1]]
  
  data(VF.trainMod8mer)
  w <- VF.trainMod8mer
  predResult <- NULL
  
  featureOut <- countSeqFeatureCpp(seqFa, w)
  featureOut_kmerCount <- featureOut$kmerCount
  seqLength <- length(seqFa)
  
  if (seqLength < 1 * 1000) {
    lasso.mod <- attr(VF.trainMod8mer, "lasso.mod_0.5k")
    rmWordID <- attr(VF.trainMod8mer, "rmWordID_0.5k")
    nullDis <- attr(VF.trainMod8mer, "nullDis_0.5k")
  }
  else if (seqLength < 3 * 1000) {
    lasso.mod <- attr(VF.trainMod8mer, "lasso.mod_1k")
    rmWordID <- attr(VF.trainMod8mer, "rmWordID_1k")
    nullDis <- attr(VF.trainMod8mer, "nullDis_1k")
  }
  else {
    lasso.mod <- attr(VF.trainMod8mer, "lasso.mod_3k")
    rmWordID <- attr(VF.trainMod8mer, "rmWordID_3k")
    nullDis <- attr(VF.trainMod8mer, "nullDis_3k")
  }
  
  lasso.pred <- predict(lasso.mod, t(as.matrix(featureOut_kmerCount[-rmWordID])),
                        type = "response")
  pvalue <- mean(nullDis > as.numeric(lasso.pred))
  print(paste("len", seqLength,
              "score", round(lasso.pred, 4), "pvalue", round(pvalue,
                                                             4)))
  
  predResult <- rbind(predResult, c(seqLength, lasso.pred, pvalue))
  colnames(predResult) <- c("length", "score", "pvalue")
  predResult_df <- as.data.frame(predResult)
  predResult_df$length <- as.numeric(as.character(predResult_df$length))
  predResult_df$score <- as.numeric(as.character(predResult_df$score))
  predResult_df$pvalue <- as.numeric(as.character(predResult_df$pvalue))
  predResult_df
}
environment(parVF.run) <- asNamespace('VirFinder')

parVF.pred <- function(inFaFile, taskcpus) {
  dnastringset <- readDNAStringSet(inFaFile)
  predResult_df <- do.call(rbind, mclapply(dnastringset, parVF.run, mc.preschedule=F, mc.cores=taskcpus))
  predResult_df$name <- rownames(predResult_df)
  predResult_df <- predResult_df[,c("name","length","score","pvalue")]
  predResult_df <- predResult_df[names(dnastringset),]
  predResult_df
}
environment(parVF.pred) <- asNamespace('VirFinder')