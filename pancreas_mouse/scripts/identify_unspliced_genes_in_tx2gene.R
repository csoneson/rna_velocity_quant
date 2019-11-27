args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

print(intsv)
print(outtsv)

x <- read.delim(intsv, header = FALSE, as.is = TRUE)
intrs <- grep("-I", x$V1)
x$V2[intrs] <- paste0(x$V2[intrs], "I.")

write.table(x, file = outtsv, row.names = FALSE,
            col.names = FALSE, sep = "\t", quote = FALSE)

date()
sessionInfo()
