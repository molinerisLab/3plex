########################
#
#	ChIRP-seq RAW
#

ChIRP_raw/hg19-DRAIR_rep1.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE156121&format=file&file=GSE156121%5FTHP1%2DDRAIR%5Frep1%2D1e%2D5%5Fpeaks%2Etxt%2Egz'

ChIRP_raw/hg19-DRAIR_rep2.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE156121&format=file&file=GSE156121%5FTHP1%2DDRAIR%5Frep2%2D1e%2D5%5Fpeaks%2Etxt%2Egz'

ChIRP_raw/hg19-DACOR.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1854427&format=file&file=GSM1854427%5FDACOR1%5FChIRP%5Fanalysis%2Etxt%2Egz'

ChIRP_raw/hg19-lincDUSP.xlsx:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3073889&format=file&file=GSM3073889%5FlincDUSP%5FChIRP%2Dseq%5Fpeaks%2Exlsx'

ChIRP_raw/hg19-SLEAR.odd1_peaks.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2801160&format=file&file=GSM2801160%5Fodd1%5Fpeaks%2Etxt%2Egz'

ChIRP_raw/hg19-SLEAR.even1_peaks.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2801161&format=file&file=GSM2801161%5Feven1%5Fpeaks%2Etxt%2Egz'

ChIRP_raw/hg19-SLEAR.odd2_peaks.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2801164&format=file&file=GSM2801164%5Fodd2%5Fpeaks%2Etxt%2Egz'

ChIRP_raw/hg19-SLEAR.even2_peaks.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2801165&format=file&file=GSM2801165%5Feven2%5Fpeaks%2Etxt%2Egz'

ChIRP_raw/hg19-SLEAR.odd3_peaks.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2801168&format=file&file=GSM2801168%5Fodd3%5Fpeaks%2Etxt%2Egz'

ChIRP_raw/hg19-SLEAR.even3_peaks.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2801169&format=file&file=GSM2801169%5Feven3%5Fpeaks%2Etxt%2Egz'

ChIRP_raw/hg19-DINO.even.narrowPeak.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2011201&format=file&file=GSM2011201%5FAS1%5F1%5FDINO%5FChIRP%5FEVEN%2EnarrowPeak%2Egz'

ChIRP_raw/hg19-DINO.odd.narrowPeak.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2011202&format=file&file=GSM2011202%5FAS1%5F2%5FDINO%5FChIRP%5FODD%2EnarrowPeak%2Egz'

ChIRP_raw/hgX-SNHG11.xlsx:
	wget -O $@ 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41388-020-01512-8/MediaObjects/41388_2020_1512_MOESM2_ESM.xls'

ChIRP_raw/hg19-SRA_rep1.peaks.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4632919&format=file&file=GSM4632919%5FK562%5FSRA%5F1%5Fp0001%5Fpeaks%2Ebed%2Egz'

ChIRP_raw/hg19-SRA_rep2.peaks.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4632920&format=file&file=GSM4632920%5FK562%5FSRA%5F2%5Fp0001%5Fpeaks%2Ebed%2Egz'

ChIRP_raw/hg19-SRA_rep3.peaks.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4632921&format=file&file=GSM4632921%5FK562%5FSRA%5F3%5Fp0001%5Fpeaks%2Ebed%2Egz'

ChIRP_raw/hg19-CDKN2B-AS1.peaks.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4947063&format=file&file=GSM4947063%5FANRIL%5FHEK293%5Fpeaks%5Fhighscore15%5Frep1%2Ebed%2Egz'

ChIRP_raw/hg38-NEAT1_TFR1.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3417018&format=file&file=GSM3417018%5FCaptureSeq%5FNEAT1%5FTFR1%5Fpeaks%2Ebed%2Egz'

ChIRP_raw/hg19-MEG3.xlsx:
	wget -O $@ 'https://static-content.springer.com/esm/art%3A10.1038%2Fncomms8743/MediaObjects/41467_2015_BFncomms8743_MOESM364_ESM.xls'

ChIRP_raw/mm9-Bloodlinc.bed.gz:
	wget -O $@ 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97119/suppl/GSE97119_EC7_consensus_aboveThresholds_peaks.bed.gz'

ChIRP_raw/mm10-Gm4665.txt.gz:
	wget -O $@ 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145189/suppl/GSE145189%5FMACS2%5FLncRNA%5FInput%2Etxt%2Egz'

ChIRP_raw/mm10-1200007C13Rik@4T1.xlsx:
	wget -O $@ 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-20207-y/MediaObjects/41467_2020_20207_MOESM5_ESM.xlsx'

ChIRP_raw/mm9-9530018H14Rik@GLIAL.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70986&format=file&file=GSE70986%5Fglial%2Dpeaks%2Ebed%2Egz'

ChIRP_raw/mm9-9530018H14Rik@LIMB.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE70986&format=file&file=GSE70986%5FLimb%2Dpeaks%2Ebed%2Egz'

ChIRP_raw/mm10-AK156552.1.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE73805&format=file&file=GSE73805%5FPanct1%2Ebed%2Egz'

ChIRP_raw/mm10-1200007C13Rik@NDL.xlsx:
	wget -O $@ 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-20207-y/MediaObjects/41467_2020_20207_MOESM6_ESM.xlsx'

ChIRP_raw/mm10-Tug1.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE77493&format=file&file=GSE77493%5Fmasterlist%2Ebed%2Etxt%2Egz'

#########################
#
#	ChIRP-seq
#


ChIRP_data/hg18-HOTAIR.bed.gz:
	wget -O $@ 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31332/suppl/GSE31332%5Fhotair%5Foe%5Fpeaks%2Ebed%2Egz'

ChIRP_data/hg18-TERC.bed.gz:
	wget -O $@ 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31332/suppl/GSE31332%5Fterc%5Fpeaks%2Ebed%2Egz'

ChIRP_data/hg19-HOTTIP.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3161929&format=file&file=GSM3161929%5FWT%2DCHIRP%2Dseq%5Fpeaks%2Ebed%2Egz'

ChIRP_data/hg19-LincASEN.bed.gz:	
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3673641&format=file&file=GSM3673641%5F713%5F2%5Fp01%5Fpeaks%2Ebed%2Egz'

ChIRP_data/hg19-SRA.bed.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1415922&format=file&file=GSM1415922%5FNTERA2%5FSRA%2Ebed%2Egz'

ChIRP_raw/mm9-Rn7sk.xlsx:
	wget -O $@ 'https://static-content.springer.com/esm/art%3A10.1038%2Fnsmb.3176/MediaObjects/41594_2016_BFnsmb3176_MOESM31_ESM.xlsx'




############################
#
#	Raw peaks processing
#


ChIRP_narrowPeaks/hgX-SNHG11.odd.narrowPeak.gz: ChIRP_raw/hgX-SNHG11.xlsx
	xls2tab $< | fasta2tab | bawk '$$1=="ODD probes"' | cut -f2- | unhead | sort -k1,1 -k2,2n | uniq | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,"NA",".",$$8,-1,$$9,$$5}' | gzip > $@
ChIRP_narrowPeaks/hgX-SNHG11.even.narrowPeak.gz: ChIRP_raw/hgX-SNHG11.xlsx
	xls2tab $< | fasta2tab | bawk '$$1=="EVEN Probes"' | cut -f2- | unhead | sort -k1,1 -k2,2n | uniq | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,"NA",".",$$8,-1,$$9,$$5}' | gzip > $@
ChIRP_narrowPeaks/hg19-DACOR.narrowPeak.gz: ChIRP_raw/hg19-DACOR.txt.gz
	zcat $< | unhead | sort -k1,1 -k2,2n | cut -f1,2,3,5,14,15,19,20 | uniq | bawk '{print "chr"$$1,$$2,$$3,"chirp_peak_"NR,$$14,$$5,$$15,$$17,$$19,$$20}' | gzip > $@
ChIRP_narrowPeaks/hg19-DACOR.narrowPeak.gz: ChIRP_raw/hg19-DACOR.txt.gz
	zcat $< | unhead | sort -k1,1 -k2,2n | bawk '{print "chr"$$1,$$2,$$3,"chirp_peak_"NR,$$14,$$5,$$15,$$19,$$20}' | gzip > $@
ChIRP_narrowPeaks/hg19-DRAIR_%.narrowPeak.gz: ChIRP_raw/hg19-DRAIR_%.txt.gz
	zgrep -v '^#' $< | unhead -n 2 | sort -k1,1 -k2,2n | uniq | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,"NA",".",$$8,$$7,$$9,$$5-$$2}' | gzip > $@
ChIRP_narrowPeaks/hg19-DINO%.narrowPeak.gz: ChIRP_raw/hg19-DINO%.narrowPeak.gz
	cp $< $@
ChIRP_narrowPeaks/hg19-SLEAR.%.narrowPeak.gz: ChIRP_raw/hg19-SLEAR.%_peaks.txt.gz
	zcat $< | unhead | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,"NA",".",$$8,$$7,$$9,$$5-$$2}' | gzip > $@
ChIRP_narrowPeaks/hg19-CDKN2B-AS1.narrowPeak.gz: ChIRP_raw/hg19-CDKN2B-AS1.peaks.bed.gz
	zcat $< | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,$$5~9,-1}' | gzip > $@
ChIRP_narrowPeaks/hg19-SRA_%.narrowPeak.gz: ChIRP_raw/hg19-SRA_%.peaks.bed.gz
	zcat $< | dos2unix | unhead | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,"NA",".","NA",$$5,-1,"NA"}' | gzip > $@
ChIRP_narrowPeaks/hg38-NEAT1_TFR1.narrowPeak.gz: ChIRP_raw/hg38-NEAT1_TFR1.bed.gz
	zcat $< | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,"NA",$$6,$$5,-1,-1,-1}' | gzip > $@

.META: %.narrowPeaks.gz
	1	chr	chr1
	2	start	5790513
	3	end	5790724
	4	peak_name	chirp_peak_1
	5	score	NA
	6	strand	.
	7	logFC	5.90736
	8	-log10pvalue	
	9	-log10qvalue	4.60985
	10	relative_summit	79


GSM3073889_lincDUSP_ChIRP-seq_peaks-2.csv:
	@echo 'Uploaded from local file'
ChIRP_narrowPeaks/hg19-lincDUSP.narrowPeak.gz: GSM3073889_lincDUSP_ChIRP-seq_peaks-2.csv
	tr ";" "\t" < $< | dos2unix | unhead | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,"NA",".",$$6,$$5,-1,"NA"}' | gzip > $@

ChIRP_narrowPeaks/%.narrowPeak.gz: ChIRP_data/%.bed.gz
	zcat $< | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,"NA",".","NA",-1,-1,"NA"}' | gzip > $@

ChIRP_narrowPeaks/hg19-MEG3.narrowPeak.gz: ChIRP_raw/hg19-MEG3.xlsx
	xls2tab $< | unhead -n 2 | cut -f2-4 | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR,"NA",".","NA",-1,-1,"NA"}' | gzip > $@

ChIRP_narrowPeaks/mm%.narrowPeaks.gz: ChIRP_raw/mm%.bed.gz
	zcat $< | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR}' | gzip > $@

ChIRP_narrowPeaks/mm%.narrowPeaks.gz: ChIRP_raw/mm%.xlsx
	xls2tab $< | unhead | bawk '{print $$2,$$3,$$4,"chirp_peak_"NR}' | gzip > $@



#############################
#
#	Converting between genome versions
#

ChIRP_narrowPeaks/hg38-%.bed: ChIRP_narrowPeaks/hg19-%.narrowPeak.gz /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/hg19ToHg38.over.chain.gz
	liftOver <(zcat $< | cut -f1-4) $^2 $@ $@.unmap

ChIRP_narrowPeaks/hg38-%.bed: ChIRP_narrowPeaks/hg18-%.narrowPeak.gz /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/hg18ToHg38.over.chain.gz
	liftOver <(zcat $< | cut -f1-4) $^2 $@ $@.unmap

ChIRP_narrowPeaks/hg38-NEAT1.bed: ChIRP_narrowPeaks/hg38-NEAT1_TFR1.narrowPeak.gz
	zcat $< | cut -f1-4 > $@

ChIRP_narrowPeaks/mm10-%.bed: ChIRP_narrowPeaks/mm9-%.narrowPeaks.gz /sto1/ref/bioinfotree/task/gencode/dataset/mmusculus/M25/mm9ToMm10.over.chain.gz
	liftOver <(zcat $< | cut -f1-4) $^2 $@ $@.unmap



###############################
#
#	Peaks intersection
#

ChIRP_narrowPeaks/hg38-SLEAR.intersection6.bed: ChIRP_narrowPeaks/hg38-SLEAR.even1.bed ChIRP_narrowPeaks/hg38-SLEAR.even2.bed ChIRP_narrowPeaks/hg38-SLEAR.even3.bed ChIRP_narrowPeaks/hg38-SLEAR.odd1.bed ChIRP_narrowPeaks/hg38-SLEAR.odd2.bed ChIRP_narrowPeaks/hg38-SLEAR.odd3.bed
	matrix_reduce -t 'ChIRP_narrowPeaks/hg38-SLEAR.*.bed' -l '$^' | bawk '{print $$2,$$3,$$4,$$1}' | bedtools sort | union --allow-duplicates | expandsets 4 | collapsesets 4 | tr ";" "\t" | bawk 'NF==9 {print $$1~3,"chirp_peak_" NR}' > $@

ChIRP_narrowPeaks/hg38-SRA.intersection4.bed: ChIRP_narrowPeaks/hg38-SRA.bed ChIRP_narrowPeaks/hg38-SRA_rep1.bed ChIRP_narrowPeaks/hg38-SRA_rep2.bed ChIRP_narrowPeaks/hg38-SRA_rep3.bed
	matrix_reduce -t 'ChIRP_narrowPeaks/hg38-SRA*' -l '$^' | bawk '{print $$2,$$3,$$4,$$1}' | bedtools sort | union --allow-duplicates | expandsets 4 | collapsesets 4 | tr ";" "\t" | bawk 'NF==6 {print $$1~3,"chirp_peak_" NR}' > $@

ChIRP_narrowPeaks/hg38-DINO.intersection2.bed: ChIRP_narrowPeaks/hg38-DINO.even.bed ChIRP_narrowPeaks/hg38-DINO.odd.bed
	matrix_reduce -t 'ChIRP_narrowPeaks/hg38-DINO.*' | bawk '{print $$2,$$3,$$4,$$1}' | bedtools sort | union --allow-duplicates | expandsets 4 | collapsesets 4 | tr ";" "\t" | bawk 'NF==5 {print $$1~3,"chirp_peak_" NR}' > $@

ChIRP_narrowPeaks/hg38-DRAIR.intersection2.bed: ChIRP_narrowPeaks/hg38-DRAIR_rep1.bed ChIRP_narrowPeaks/hg38-DRAIR_rep2.bed
	matrix_reduce -t 'ChIRP_narrowPeaks/hg38-DRAIR*' | bawk '{print $$2,$$3,$$4,$$1}' | bedtools sort | union --allow-duplicates | expandsets 4 | collapsesets 4 | tr ";" "\t" | bawk 'NF==5 {print $$1~3,"chirp_peak_" NR}' > $@





##############################
#
#	Total peaks
#

ChIRP_narrowPeaks.bed: ChIRP_narrowPeaks/hg38-MEG3.bed ChIRP_narrowPeaks/hg38-CDKN2B-AS1.bed ChIRP_narrowPeaks/hg38-NEAT1.bed ChIRP_narrowPeaks/hg38-DACOR.bed ChIRP_narrowPeaks/hg38-DINO.intersection2.bed ChIRP_narrowPeaks/hg38-DRAIR.intersection2.bed ChIRP_narrowPeaks/hg38-HOTAIR.bed ChIRP_narrowPeaks/hg38-HOTTIP.bed ChIRP_narrowPeaks/hg38-lincDUSP.bed ChIRP_narrowPeaks/hg38-SLEAR.intersection6.bed ChIRP_narrowPeaks/hg38-SRA.intersection4.bed ChIRP_narrowPeaks/hg38-TERC.bed
	matrix_reduce -t 'ChIRP_narrowPeaks/hg38-*.bed' -l '$^' | bawk '{print $$2~5";"$$1}' \
        | perl -pe 's/lincDUSP/LINC01605/; s/DRAIR/CPEB2-DT/; s/DACOR/AC087482.1/; s/DINO/Z85996.1/; s/.intersection\d+//; s/SRA/SRA1/; s/SLEAR/AC018781.1/' > $@

ChIRP.bed.fa: /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/GRCh38.primary_assembly.genome.fa ChIRP_narrowPeaks.bed
	bedtools getfasta -name -fi $< -bed $^2 > $@


##############################
#
#	lncRNA selection
#

v14_selected_ssRNA: /sto1/ref/bioinfotree/task/encode/dataset/20210615/standard_processing/dataset/v1/top_lncRNA.ALLconditions
	cp -a $< $@

cCRE-%.bed: /sto1/ref/bioinfotree/task/encode-screen/dataset/v13/hg38-%.bed
	cp -a $< $@



# ---------------- #
#  ReChIRP project #
# ---------------- #

prj_lncRNA.dict:
	@echo 'https://docs.google.com/spreadsheets/d/13o1DAEC3XndtZ3ACe-3QvQXECJIgh0de/edit#gid=348805968'

### IDR HUMAN ###

ReChIRP/idr/%.idr.conservative_peak.regionPeak.gz: ReChIRP/%.idr.conservative_peak.regionPeak.gz
	mkdir -p `dirname $@`; \
	bawk '{print $$1~3,"chirp_peak_" NR}' $< | gzip > $@

ReChIRP/idr/idr.conservative_peak.regionPeak.human_matrix.gz:
	matrix_reduce -t 'ReChIRP/idr/*.idr.conservative_peak.regionPeak.gz' | bawk '{print $$2~5,$$1}' | gzip > $@


### OVERLAP HUMAN ###

ReChIRP/overlap/AC018781.9.overlap.regionPeak.qTh005.sorted.gz: /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/GSE104479_PRJNA412810_rep1/overlap.regionPeak.qTh005.sorted.gz /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/GSE104479_PRJNA412810_rep2/overlap.regionPeak.qTh005.sorted.gz /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/GSE104479_PRJNA412810_rep3/overlap.regionPeak.qTh005.sorted.gz
	multiIntersectBed -i $< $^2 $^3 | bawk '$$4==3 {print $$1~3,"chirp_peak_" NR}' | gzip > $@

ReChIRP/overlap/%.overlap.regionPeak.qTh005.sorted.gz: ReChIRP/%.overlap.regionPeak.qTh005.sorted.gz
	mkdir -p `dirname $@`; \
        bawk '{print $$1~3,"chirp_peak_" NR}' $< | gzip > $@

ReChIRP/overlap/overlap.regionPeak.qTh005.sorted.human_matrix.gz:
	matrix_reduce -t 'ReChIRP/overlap/*.overlap.regionPeak.qTh005.sorted.gz' | bawk '{print $$2~5,$$1}' | gzip > $@

ReChIRP/overlap/overlap.conservative_peak.regionPeak.qTh005.human_matrix.gz: GSE_ssRNA_overlap /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/GSE104479_PRJNA412810_rep1/overlap.conservative_peak.regionPeak.rep_123.qTh005.gz
	matrix_reduce -t '/sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/*/peak/overlap_reproducibility/overlap.conservative_peak.regionPeak.gz' | bawk '$$10>-log(0.05)/log(10)' | translate -k $< 1 | bawk '{print $$2~4,"chirp_peak_" NR,$$1}' > $@.tmp; \
	cat <(bawk '{print $$1~3,"chirp_peak_" NR,"AC018781.1"}' $^2) $@.tmp | gzip > $@; \
	rm $@.tmp



### IDR + OVERLAP TOP 1000 HUMAN ###

ReChIRP/idr_overlap_top1000/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.top1000.gz: GSE_ssRNA_idr-overlap /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.top1000.gz /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/GSE104479_PRJNA412810_rep1/idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.top1000.rep_123.gz
	mkdir -p `dirname $@`; \
	zcat $^2 | translate -k $< 1 | bawk '{print $$2~4,"chirp_peak_" NR,$$1}' > $@.tmp; \
	cat <(bawk '{print $$1~3,"chirp_peak_" NR,"AC018781.1"}' $^3) $@.tmp | gzip > $@; \
	rm $@.tmp

ReChIRP/idr_overlap_top1000/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.gz: /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.gz GSE_ssRNA_idr-overlap
	zcat $< | translate -k $^2 1 | gzip > $@



### IDR MOUSE ###

mmusculus/ssRNA_fasta/%.fa:
	wget -O $@.tmp "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$*&rettype=fasta"; \
	bawk '{if(NR==1){split($$1,a," "); print a[1]} else {print $$0}}' $@.tmp > $@; \
	rm $@.tmp;
	@echo 'check if last line of the file is empty'

ReChIRP/mmusculus/idr/idr.conservative_peak.regionPeak.mouse_matrix.gz: GSE_ssRNA_mouse_idr
	mkdir -p `dirname $@`; \
	matrix_reduce -t '/sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/*/peak/idr_reproducibility/idr.conservative_peak.regionPeak.gz' | translate -k $< 1 | bawk '{print $$2~4,"chirp_peak_" NR,$$1}' | gzip > $@


### OVERLAP MOUSE ###

ReChIRP/mmusculus/overlap/overlap.conservative_peak.regionPeak.qTh005.mouse_matrix.gz: GSE_ssRNA_mouse_overlap
	mkdir -p `dirname $@`; \
	matrix_reduce -t '/sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/*/peak/overlap_reproducibility/overlap.conservative_peak.regionPeak.gz' | bawk '$$10>-log(0.05)/log(10)' | translate -k $< 1 | bawk '{print $$2~4,"chirp_peak_" NR,$$1}' | gzip > $@


### IDR + OVERLAP TOP 1000 MOUSE ###

ReChIRP/mmusculus/idr_overlap_top1000/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.top1000.gz: GSE_ssRNA_mouse_idr-overlap /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.top1000.gz
	mkdir -p `dirname $@`; \
        zcat $^2 | translate -k $< 1 | bawk '{print $$4~6,"chirp_peak_" NR,$$2}' | gzip > $@

ReChIRP/mmusculus/idr_overlap_top1000/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.gz: /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.gz GSE_ssRNA_mouse_idr-overlap
	zcat $< | translate -k $^2 1 | gzip > $@


### IDR OPTIMAL HUMAN ###

ReChIRP/idr/idr.optimal_peak.regionPeak.human_matrix.gz: GSE_ssRNA_idr_optimal /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/GSE104479_PRJNA412810_rep1/idr.optimal_peak.regionPeak.rep_123.gz
	mkdir -p `dirname $@`; \
        matrix_reduce -t '/sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/*/peak/idr_reproducibility/idr.optimal_peak.regionPeak.gz' | translate -k $< 1 | bawk '{print $$2~4,"chirp_peak_" NR,$$1}' > $@.tmp; \
	cat <(bawk '{print $$1~3,"chirp_peak_" NR,"AC018781.1"}' $^2) $@.tmp | gzip > $@; \
        rm $@.tmp


### IDR OPTIMAL MOUSE ###

ReChIRP/mmusculus/idr/idr.optimal_peak.regionPeak.mouse_matrix.gz: GSE_ssRNA_mouse_idr_optimal
	mkdir -p `dirname $@`; \
	matrix_reduce -t '/sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/*/peak/idr_reproducibility/idr.optimal_peak.regionPeak.gz' | translate -k $< 1 | bawk '{print $$2~4,"chirp_peak_" NR,$$1}' | gzip > $@


### ALL POS-NEG PEAKS ###
ALL_v8.neg_pos.bed.gz:
	matrix_reduce -t '../../../dataset/*/*.neg_pos_rand.bed' | grep '^v8' | tr ";" "\t" | gzip > $@


CLEAN_all_reproducibility.regionPeak.counts_matrix.overlap_custom.overlap_conservative_qTh005.tsv:
	@echo 'Downloaded from https://docs.google.com/spreadsheets/d/1ZOQSaY3k7ucJFEFf633jABl8KLdIsdRp/edit#gid=348805968'

CLEAN_all_reproducibility.regionPeak.counts_matrix.overlap_custom.overlap_conservative_qTh005.bed_peaks.tsv: /sto1/epigen/TPXcCRE/dataset/v8_ChIRP_neg_rand/ALL_neg_pos_rand.bed CLEAN_all_reproducibility.regionPeak.counts_matrix.overlap_custom.overlap_conservative_qTh005.tsv
	translate -a -r -v -e NA <(cat $< | tr ";" "\t" | cut -f5 | symbol_count) 2 < $^2 > $@

ChIRP_narrowPeaks/mm%.narrowPeaks.gz: ChIRP_raw/mm%.bed.gz
	zcat $< | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR}' | gzip > $@

ChIRP_narrowPeaks/mm10-Gm4665.narrowPeaks.gz: ChIRP_raw/mm10-Gm4665.txt.gz
	zcat $< | unhead | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR}' | gzip > $^2
ChIRP_narrowPeaks/mm10-AK156552.1.narrowPeaks.gz: ChIRP_raw/mm10-AK156552.1.bed.gz
	zcat $< | sort -k1,1 -k2,2n | bawk '{print $$1,$$2,$$3,"chirp_peak_"NR}' | gzip > $^2

ChIRP_narrowPeaks/mm%.narrowPeaks.gz: ChIRP_raw/mm%.xlsx
	xls2tab $< | unhead | bawk '{print $$2,$$3,$$4,"chirp_peak_"NR}' | gzip > $@
mm10.ChIRP_narrowPeaks.gz:
	matrix_reduce -t 'ChIRP_narrowPeaks/mm10-*.narrowPeaks.gz' | bawk '{print $$2~5";"$$1}' | gzip > $@
mmusculus/ssRNA_fasta/lncSmad7_parts.fa: /sto1/epigen/Maldotti_lncSMAD7/lncSMAD7_target/local/share/data/mouse_lncRNASmad7.fa
	cp $< $@

all_reproducibility.all_bed.regionPeak.counts_matrix: all_bed.regionPeak.counts_matrix all_reproducibility.regionPeak.counts_matrix
	translate -a -r -v -e NA $< 2 < $^2 > $@

all_bed.regionPeak.counts_matrix: /sto1/epigen/TPXcCRE/dataset/v8_ChIRP_neg_rand/ChIRP.bed.split mm10.ChIRP_narrowPeaks.gz
	cut -f5 $< | symbol_count | bawk 'BEGIN{print "lncRNA","available_bed"} {print}' > $@; \
	zcat $< | tr ";" "\t" | cut -f5 | symbol_count >> $@



### RADICL-seq ###

RADICLseq/GSE132190_mESC_NPM_n1_significant.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132190&format=file&file=GSE132190%5FmESC%5FNPM%5Fn1%5Fsignificant%2Etxt%2Egz'

RADICLseq/GSE132190_mESC_NPM_n2_significant.txt.gz:
	wget -O $@ 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132190&format=file&file=GSE132190%5FmESC%5FNPM%5Fn2%5Fsignificant%2Etxt%2Egz'

#RADICLseq/GSE132190_mESC_NPM_n%_significant.merge_distinct.gz: RADICLseq/GSE132190_mESC_NPM_n%_significant.txt.gz /home/reference_data/bioinfotree/task/gencode/dataset/mmusculus/M25/chrom.info
#	bawk '$$29!="."' $< | bedtools slop -i - -g $^2 -b 2500 | bedtools sort | bedtools merge -i - -c 29 -o distinct | gzip > $@

#RADICLseq/GSE132190_mESC_NPM_ALL_significant.merge_distinct.matrix_reduce.merge_replicates.biotypes.gz: /home/reference_data/bioinfotree/task/gencode/dataset/mmusculus/M25/primary_assembly.annotation.ensg2gene_symbol2biotype.map
#	matrix_reduce -t 'RADICLseq/GSE132190_mESC_NPM_*_significant.merge_distinct.gz' | \
#	expandsets 5 -s , | \
#	bawk '{print $$2~5";"$$1}' | \
#	bedtools sort | bedtools merge -i - -c 4 -o distinct | \
#	expandsets 4 -s , | tr ";" "\t" | collapsesets 5 | \
#	translate -a -k -d <(cut -f2,3 $<) 4 | \
#	gzip > $@

.META: GSE132190_mESC_NPM_n*_significant.txt.gz
	1	dna_chr	chr3
	2	dna_b	108117864
	3	dna_e	108117864
	4	dna_name	.
	5	dna_score	.
	6	dna_strand	-
	7	dna_geneID	ENSMUSG00000000001.4
	8	dna_gene_biotype_1	protein_coding
	9	dna_gene_biotype_2	protein_coding
	10	dna_gene_feature	intron
	11	rna_chr	chr3
	12	rna_b	108117866
	13	rna_e	108117866
	14	rna_name	chr3_4325
	15	rna_x	1
	16	controlID	RADICL_mESC_Control
	17	inter_pval	1.0383050834521e-60
	18	inter_qval	3.13991552556332e-57
	19	repeat	.
	20	repeat	-1
	21	repeat	-1
	22	repeat	.
	23	repeat	.
	24	repeat	.
	25	repeat	.
	26	rna_gene_chr	chr3
	27	rna_gene_b	108107279
	28	rna_gene_e	108146146
	29	rna_gene_name	Gnai3
	30	rna_geneID	ENSMUSG00000000001.4
	31	rna_gene_strand	-
	32	rna_gene_biotype	protein-coding

