args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
  library(Biostrings)
})

print(infa)
print(outfa)

x <- readDNAStringSet(infa)
idx1 <- grep("-I", names(x))
idx2 <- grep("-I", names(x), invert = TRUE)
x <- x[c(idx1, idx2)]

writeXStringSet(x, file = outfa, compress = TRUE)

date()
sessionInfo()
