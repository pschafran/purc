print(format(Sys.time()))

# load libraries
library(dada2);packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(gridExtra); packageVersion("gridExtra")

path <- "."

fns <- list.files(path, pattern="fastq", full.names=TRUE)
Fprimer <- "GATCTTTATGAACAATGTGGGA"
Rprimer <- "GAAATACCTGATTTGTAACCTA"
rc <- dada2:::rc
print("Forward Primer:")
print(Fprimer)
print("Reverse Primer:")
print(Rprimer)

nops <- file.path(path, "noprimers", basename(fns))
prim <- removePrimers(fns, nops, primer.fwd=Fprimer, primer.rev=dada2:::rc(Rprimer), orient=TRUE)
lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)

minLen <- 0
maxLen <- 0
maxEE <- 10

# Calculate outlier boundaries if min/max length not specified in config file
if (minLen == 0){
iqr <- IQR(lens)
q25 <- quantile(lens, 0.25)
minLen <- q25-(1.5*iqr)
}
if (maxLen == 0){
iqr <- IQR(lens)
q75 <- quantile(lens, 0.75)
maxLen <- q75+(1.5*iqr)
}

print("Minimum length:")
print(minLen)
print("Maximum length:")
print(maxLen)
print("Maximum expected errors:")
print(maxEE)

pdf("Isoetes_straightLeaves_Taylor6989-1_2_read_lengths.pdf")
hist(lens, 100, xlim = c(min(lens)-1, max(lens)+1), xlab = "Length (bp)", main = "Isoetes_straightLeaves_Taylor6989-1_2")
abline(v= c(minLen, maxLen), lty=c(2,2))
dev.off()

filts <- file.path(path, "noprimers", "filtered", basename(fns))
track <- filterAndTrim(nops, filts, minQ=3, minLen=minLen, maxLen=maxLen, maxN=0, rm.phix=FALSE, maxEE=maxEE)
print("Reads filtered:")
print(track)

drp <- derepFastq(filts, verbose=TRUE)
err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
pdf("Isoetes_straightLeaves_Taylor6989-1_2_error_profile.pdf")
plotErrors(err)
dev.off()

dd2 <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)
st <- makeSequenceTable(dd2)

# Identify chimeras
st.nochim <- removeBimeraDenovo(st, method="consensus", multithread=TRUE, verbose=TRUE)
print("Non-chimeric reads:")
print(dim(st.nochim))
print("Percent non-chimeric reads:")
print(sum(st.nochim)/sum(st)*100)

print("SUMMARY")
print(cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2], denoised=sum(dd2$denoised), ASVs=dim(st.nochim)[2]))

# Output final ASV seqs
outputNames <- vector()
for (i in 1:dim(st.nochim)[2]){
  newName <- paste("Isoetes_straightLeaves_Taylor6989-1_2_ASV", i, sep = "")
  outputNames <- c(outputNames, newName)
}
uniquesToFasta(st.nochim, "Isoetes_straightLeaves_Taylor6989-1_2_ASVs.fasta", ids=outputNames)
print(format(Sys.time()))
quit()
