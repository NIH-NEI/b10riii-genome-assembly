# R script for Jellyfish output analysis
# Author : Vijay Nagarajan PhD
# NIH/NEI/LI

# read jelly histo
jel<-read.table("jellyhisto.txt")

# plot histogram of kmers
pdf("jel.pdf")
plot(jel[3:100,],type="l")
dev.off()

# plot histogram with points
pdf("jel2.pdf")
plot(jel[3:100,],type="l")
points(jel[5:100,])
dev.off()

# calculate total kmers
sum(as.numeric(jel[5:9989,1]*jel[5:9989,2]))

# look at around peak and identify the exact peak - 18
jel[10:30,]

# calculate the genome size
sum(as.numeric(jel[5:9989,1]*jel[5:9989,2]))/18

# plot poisson
pdf("jelpoi.pdf")
singleC <- sum(as.numeric(jel[5:30,1]*jel[5:30,2]))/18
plot(3:100,dpois(3:100, 18)*singleC, type = "l", col=3, lty=2)
lines(jel[3:100,],type="l")
dev.off()

