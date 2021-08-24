################
### function ###
################

# library(multtest)
library(GO.db)


# ng.BHFDR ------------
#p valueのベクトルから多重検定の補正をする関数のラッパー 
#Benjamini & Hochberg (1995) step-up FDR
#FDR補正後のp-valueを入力したp-valueに対応した順で返す 
# 170821： multtestが無くなった？ので、p.adjustに置き換え

#ng.BHFDR <- function(rawp){
#  tmp <- mt.rawp2adjp(rawp, proc=c("BH"))
#  adjp <- tmp$adjp
#  index <- tmp$index
#  out <- adjp[order(index), 2]
#  return(out)
#}


# multiple fisher test --------------------------
ng.mft <- function(
  cgt, #output from ng.MakeContigGOidTable, [,1]:"locus", [,2]:"GOid"
  gn.test, #contig names for test
  alternative="greater"
){
  
  #cat(sprintf("%s\n", Sys.time()))
  
  gid.u <- unique(cgt[,"GOid"])
  
  ft.in <- matrix(0, nrow=length(gid.u), ncol=9)
  colnames(ft.in) <- c("xtt", "xft", "xtf", "xff", "xnt", "xnf", "xtn", "xfn", "xnn")
  rownames(ft.in) <- gid.u
  
  #               gn.test
  #             TRUE FALSE
  #Group  TRUE   xtt   xft   xnt
  #      FALSE   xtf   xff   xnf
  #              xtn   xfn   xnn
  
  ft.in[,"xnn"] <- length(unique(cgt[, "locus"]))
  
  gn.pp.gid <- table(cgt[, "GOid"])
  ft.in[names(gn.pp.gid), "xnt"] <- gn.pp.gid
  ft.in[,"xnf"] <- ft.in[,"xnn"] - ft.in[,"xnt"]
  
  ft.in[,"xtn"] <- length(intersect(gn.test, unique(cgt[, "locus"])))
  ft.in[,"xfn"] <- ft.in[,"xnn"] - ft.in[,"xtn"]
  
  gsea.test <- cgt[is.element(cgt[,"locus"], gn.test), ]
  gn.test.gid <- table(gsea.test[, "GOid"])
  ft.in[names(gn.test.gid), "xtt"] <- gn.test.gid
  
  ft.in[,"xtf"] <- ft.in[,"xtn"] - ft.in[,"xtt"]
  ft.in[,"xft"] <- ft.in[,"xnt"] - ft.in[,"xtt"]
  ft.in[,"xff"] <- ft.in[,"xnf"] - ft.in[,"xtf"]
  
  #cat(sprintf("%s\n", Sys.time()))
  
  #Fisher's exact test.  8? sec
  fr <- rep(1, nrow(ft.in))
  dt <- rep(1, nrow(ft.in))
  for(i in 1:nrow(ft.in)){
    start <- Sys.time()
    if(ft.in[i,"xtn"] > 1 && ft.in[i,"xnt"] > 1){ 
      contable <- matrix(ft.in[i, 1:4], ncol=2)
      tmp <- fisher.test(contable, alternative = alternative)
      fr[i] <- tmp$p.value
    } else {
    }
    end <- Sys.time()
    dt[i] <- end - start
  }
  
  out <- cbind(fr, ft.in, dt)
  colnames(out) <- c("p.value", colnames(ft.in), "time")
  rownames(out) <- rownames(ft.in)
  
  #cat(sprintf("%s\n", Sys.time()))
  
  return(out)
  
}


# get GO terms -----------------------------
ng.GetGOTerms <- function(GOid){
  
  out <- NULL
  for(i in GOid){
    tmp <- try(get(i, GOTERM), silent=TRUE)
    if(class(tmp)=="try-error"){
      out <- c(out, "NA")
    } else {
      out <- c(out, Term(tmp))
    }
  }
  return(out)
}


# prep. GO test output table ------------------------------
ng.prepGOtestOutTable <- function(r, alpha){
  
  adp <- p.adjust(r[,"p.value"])
  
  tmp.id <- rownames(r)[adp < alpha]
  tmp.adp <- adp[adp < alpha]
  tmp.description <- ng.GetGOTerms(tmp.id)
  tmp.xnn <- r[adp < alpha, c("xtt", "xtn", "xnt", "xnn")]
  
  out <- cbind(tmp.adp, tmp.id, tmp.description, tmp.xnn)
  colnames(out) <- c("Adjusted P value", "ID", "Description",  "A & B", "A",	"B", "U")
  
  return(out)
}
