# Stephen Turner
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Last updated: Tuesday, April19, 2011
# R code for making manhattan plots and QQ plots from plink output files. 
# manhattan() with GWAS data this can take a lot of memory, recommended for use on 64bit machines only, for now. 
# Altnernatively, use bmanhattan() , i.e., base manhattan. uses base graphics. way faster.


## This is for testing purposes.
# set.seed(42)
# nchr=23
# nsnps=1000
# d=data.frame(
#     SNP=sapply(1:(nchr*nsnps), function(x) paste("rs",x,sep='')),
#     CHR=rep(1:nchr,each=nsnps), 
#     BP=rep(1:nsnps,nchr), 
# 	P=runif(nchr*nsnps)
# )
# annotatesnps <- d$SNP[7550:7750]

# manhattan plot using base graphics
manhattan <- function(dataframe, colors=c("green", "blue", "purple"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(.0001), genomewideline=-log10(5e-8), annotate=NULL, ...) {

    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    
    if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
    if (ymax=="max") ymax<-ceiling(max(d$logp))
    if (ymax<7) ymax<-7
    
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
        for (i in unique(d$CHR)) {
          if (i==1) {
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
    		} else {
    			lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
    		}
    		ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    	}
    }
    
    if (numchroms==1) {
        with(d, plot(pos, logp, pch=20, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
    }	else {
        with(d, plot(pos, logp, pch=20, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
        axis(1, at=ticks, lab=unique(d$CHR), ...)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
            icol=icol+1
    	}
    }
    
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="black", pch=17, ...)) 
    }
    
    if (suggestiveline) abline(h=suggestiveline, col="black")
    if (genomewideline) abline(h=genomewideline, col="red")
}



qq = function(p,BH=F,CI=T,FDRthres=0.05,plot=T,...)
{
    nn = length(p)
    xx =  -log10((1:nn)/(nn+1))
    dat<-cbind(sort(p),1:nn)
    q<-(nn*dat[,1])/dat[,2] # calculate q-values from p-values
    dat<-cbind(dat,q)
    if (min(dat[,3]>FDRthres)) {
        nsnps<-0
    } else {
        nsnps<-round((sum(p<=max(dat[dat[,3]<FDRthres,1]))/nn)*100,4)
    }
    if (plot) {
        plot( xx,  -sort(log10(p)),
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        cex.lab=1,mgp=c(2,1,0),pch=20,
        ... )
        if(CI)
        {
            ## create the confidence intervals
            c95 <- rep(0,nn)
            c05 <- rep(0,nn)
            ## the jth order statistic from a
            ## uniform(0,1) sample
            ## has a beta(j,n-j+1) distribution
            ## (Casella & Berger, 2002,
            ## 2nd edition, pg 230, Duxbury)
            ## this portion was posted by anonymous on
            ## http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
    
            for(i in 1:nn)
            {
             	c95[i] <- qbeta(0.95,i,nn-i+1)
                c05[i] <- qbeta(0.05,i,nn-i+1)
            }
            polygon(c(xx, rev(xx)), c(-log10(c95), rev(-log10(c05))),col = "grey", border = NA)
            #    lines(xx,-log10(c95),col='gray')
            #    lines(xx,-log10(c05),col='gray')
        }
        points(xx,  -sort(log10(p)),pch=20)
        #points(xx[dat[,3]<0.05],  -log10(dat[dat[,3]<0.05,1]),pch=20,col="pink")


        y<-max(-log10(p))
        #text(0,y,paste(nsnps,"% of SNPs have a q-value <= ",FDRthres,sep=""),pos=4)
    
        abline(0,1,col='black',lwd=2)
        if(BH)
        {
            abline(-log10(0.05),1, col='black',lty=2,lwd=1.5)
            abline(-log10(0.10),1, col='black',lty=3,lwd=1.5)
            abline(-log10(0.25),1, col='black',lty=4,lwd=1.5)
            legend('bottomright', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
            col=c('black','black','black'),lty=2:4, lwd=2, cex=1)
            abline(h=-log10(0.05/nn),col="red",lwd=2) ## bonferroni
        }
    } else {
        return(nsnps)
    }
}

