#! /usr/bin/env Rscript
#fig.name <- 'sty_map_bbk.pdf'
setwd('~/weilab/project/sty_ptm/fig_rra_threshold/compare_with_bbk/')
options(stringsAsFactors=FALSE)

x.threshold <- c(-3, 3)
y.threshold <- c(log10(0.05), -log10(0.05))

bbk <- read.table('BBK.txt',sep='\t',header = TRUE)
bbk <- bbk[,c(1,2,3,4)]
bbk <- bbk[!duplicated(bbk),]
colnames(bbk) <- c('gene_name','bbk_score',
                   'bbk_zlfc','bbk_rra')

lib1_d <- read.table('~/weilab/project/sty_ptm/analysis/dropout_sg_info_screen_score.txt',
                       sep='\t',header=TRUE)
lib1_d <- lib1_d[,c(6,11,12,15,17)]
lib1_d <- lib1_d[!duplicated(lib1_d),]
lib1_e <- read.table('~/weilab/project/sty_ptm/analysis/enrich_sg_info_screen_score.txt',
                     sep='\t',header=TRUE)
lib1_e <- lib1_e[,c(6,11,12,15,17)]
lib1_e <- lib1_e[!duplicated(lib1_e),]
lib1 <- rbind(lib1_d,lib1_e)
lib1$lib <- 'lib1'

lib2_d <- read.table('~/weilab/project/sty_ptm/analysis/lib2_dropout_sg_info_screen_score.txt',
                     sep='\t',header=TRUE)
lib2_d <- lib2_d[,c(4,10,12,14,15)]
lib2_d <- lib2_d[!duplicated(lib2_d),]
lib2_e <- read.table('~/weilab/project/sty_ptm/analysis/lib2_enrich_sg_info_screen_score.txt',
                     sep='\t',header=TRUE)
lib2_e <- lib2_e[,c(4,10,12,14,15)]
lib2_e <- lib2_e[!duplicated(lib2_e),]
lib2 <- rbind(lib2_d,lib2_e)
lib2$lib <- 'lib2'

sty <- rbind(lib1,lib2)
sty$log_rra <- -log10(sty$RRA)
sty[which(sty$zlfc<0),'log_rra'] <- -sty[which(sty$zlfc<0),
                                              'log_rra']
sty <- sty[!duplicated(sty),]
write.table(sty,
            '~/weilab/project/sty_ptm/analysis/sty_screen_result.txt',
            sep='\t',quote = FALSE,row.names = FALSE)


sty_bbk <- merge(sty,bbk,all.x = TRUE)
sty_bbk$bbk_log_rra <- -log10(sty_bbk$bbk_rra)
sty_bbk[which(sty_bbk$bbk_zlfc<0),'bbk_log_rra'] <- 
    -sty_bbk[which(sty_bbk$bbk_zlfc<0),'bbk_log_rra']


validate <- read.table('~/weilab/project/sty_ptm/fig/compare_with_bbk/compare_with_bbk_validate4.txt',
                       sep='\t',header=FALSE)
colnames(validate) <- c('gene_aa_id')


sty_bbk$color <- densCols(
    sty_bbk$log_rra,
    sty_bbk$bbk_log_rra,
    colramp=colorRampPalette(
        c('gray65', 'gray30', 'gray20', 'gray10')
    )
)

dD <- sty_bbk[which(sty_bbk$RRA < 0.001 &
                    sty_bbk$zlfc < 0 &
                    sty_bbk$bbk_rra < 0.05 &
                    sty_bbk$bbk_zlfc < 0),]
dN <- sty_bbk[which(sty_bbk$RRA < 0.001 &
                    sty_bbk$zlfc < 0 &
                    sty_bbk$bbk_rra >= 0.05),]
dE <- sty_bbk[which(sty_bbk$RRA < 0.001 &
                    sty_bbk$zlfc < 0 &
                    sty_bbk$bbk_rra < 0.05 &
                    sty_bbk$bbk_zlfc > 0),]

eD <- sty_bbk[which(sty_bbk$RRA < 0.001 &
                    sty_bbk$zlfc > 0 &
                    sty_bbk$bbk_rra < 0.05 &
                    sty_bbk$bbk_zlfc < 0),]
eN <- sty_bbk[which(sty_bbk$RRA < 0.001 &
                    sty_bbk$zlfc > 0 &
                    sty_bbk$bbk_rra >= 0.05),]
eE <- sty_bbk[which(sty_bbk$RRA < 0.001 &
                    sty_bbk$zlfc > 0 &
                    sty_bbk$bbk_rra < 0.05 &
                    sty_bbk$bbk_zlfc > 0),]




write.table(dD,'dD.txt',sep='\t',quote = FALSE,row.names = FALSE)
write.table(dN,'dN.txt',sep='\t',quote = FALSE,row.names = FALSE)
write.table(dE,'dE.txt',sep='\t',quote = FALSE,row.names = FALSE)
write.table(eD,'eD.txt',sep='\t',quote = FALSE,row.names = FALSE)
write.table(eN,'eN.txt',sep='\t',quote = FALSE,row.names = FALSE)
write.table(eE,'eE.txt',sep='\t',quote = FALSE,row.names = FALSE)





sty_bbk$color[sty_bbk$log_rra <= x.threshold[1]] <- '#e63946'

sty_bbk$color[sty_bbk$log_rra >= x.threshold[2]] <- '#023e8a'

x.limits <- c(-9, 12)
y.limits <- c(-24, 44)

old.mar <- par('mar')

#pdf('sty_map_bbk_nolable_validate.pdf', width=8, height=6, bg='transparent')

png("sty_map_bbk_lable_validate4.png",units="in", width=8, 
    height=6,res=2000, bg='transparent')

par(mar=c(5.1, 4.1, 4.1, 8.1), pty='s')

par(xpd=TRUE)

plot(
    c(0), c(0), pch=20,
    xlim=x.limits,
    ylim=y.limits,
#    xlab='',
#    ylab='',
#    xaxt='n',
#    yaxt='n',
    xlab='STY RRA(-log10)',
    ylab='BARBEKO RRA(-log10)',
    xaxs='i', yaxs='i'
)

## ## inconsistent
## rect(
##     xleft=c(x.limits[1], x.threshold[2]),
##     ybottom=c(y.threshold[1], y.limits[1]),
##     xright=c(x.threshold[1], x.limits[2]),
##     ytop=c(y.limits[2], y.threshold[2]),
##     col='#FCD5CE',
##     border = FALSE
## )

## ## consistent
## rect(
##     xleft=c(x.limits[1], x.threshold[2]),
##     ybottom=c(y.limits[1], y.threshold[2]),
##     xright=c(x.threshold[1], x.limits[2]),
##     ytop=c(y.threshold[1], y.limits[2]),
##     col='#D8E2DC',
##     border = FALSE
## )

## abline(h=y.threshold[1], col='red', lty=2, cex=3)
## abline(h=y.threshold[2], col='blue', lty=2, cex=3)

points(
    x=sty_bbk$log_rra,
    y=sty_bbk$bbk_log_rra,
    pch=20,
    cex=0.1,
    col=sty_bbk$color
)

sty_bbk_validate <- sty_bbk[which(sty_bbk$gene_aa_id %in% validate$gene_aa_id),]
text(x=sty_bbk_validate$log_rra,
     y=sty_bbk_validate$bbk_log_rra,
     labels=sty_bbk_validate$gene_aa_id,cex=0.2)

par(xpd=TRUE)

#legend(
#    x=x.limits[2] * 1.2, y=0,
#    legend=c('Enriched STY', 'Depleted STY'),
#    col=c('#023e8a', '#e63946'),
#    pch=20
#)

dev.off()

par(mar=old.mar)

if (Sys.info()['sysname'] == 'Linux') {
    system(
        paste(
            'convert -density 300',
            fig.name,
            gsub('pdf', 'png', fig.name)
        )
    )
}

## sessionInfo()
## R version 3.6.0 (2019-04-26)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux Server 7.3 (Maipo)

## Matrix products: default
## BLAS:   /gpfs/share/software/R/3.6.0/gcc/7.2.0/lib64/R/lib/libRblas.so
## LAPACK: /gpfs/share/software/R/3.6.0/gcc/7.2.0/lib64/R/lib/libRlapack.so

## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base

## loaded via a namespace (and not attached):
## [1] compiler_3.6.0     KernSmooth_2.23-15
####################
