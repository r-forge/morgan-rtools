library( "IBDhaploRtools" )

qibd.geno <- "tofiona/afr.afr.G.qibd"
ids.geno <- "tofiona/afr.afr.G.ids"

qibd.haplo <- "tofiona/afr.afr.H.qibd"
ids.haplo <- "tofiona/afr.afr.H.ids"

states.geno <- ibdhap.make.calls(qibd.filename = qibd.geno,
                                 ids.filename = ids.geno,
                                 cutoff = 0.8 )
states.haplo <- ibdhap.make.calls(qibd.filename = qibd.haplo,
                                  ids.filename = ids.haplo,
                                  cutoff = 0.8 )

true.geno <- ibdhap.make.true("tofiona/outnineibd.txt" )
true.haplo <- ibdhap.make.true( "tofiona/outfifteenibd.txt" )

## Site and Segment Summaries
sites.geno <- ibdhap.compare.loci( states.geno, true.geno, data.type = "g" )
segs.geno <- ibdhap.compare.segs( states.geno, true.geno, data.type = "g",
                                 seg.cutoff = 0.8)


sites.haplo <- ibdhap.compare.loci( states.haplo, true.haplo, data.type = "h" )
segs.haplo <- ibdhap.compare.segs( states.haplo, true.haplo, data.type = "h",
                                  seg.cutoff = 0.8)

### 4495 marker rows and 1000 chromosome columns
#############################################################################
