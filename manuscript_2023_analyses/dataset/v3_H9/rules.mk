REGREGION_DISTANCE?=5000
TPX_DISTANCE=500
GENCODE_DIR=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/hsapiens/32
#TRIPLEXATOR_PARAM?=-l 10 -e 20 -E 1
CORES?=32

TRIPLEXATOR_PARAM=-l 10 -L -1 -e 20 -E -1

#######################################
#
#	docker param
#

DOCKER_GROUP?=$(shell stat -c %g $(PWD))
DOCKER_DATA_DIR?=$(PRJ_ROOT)

ssRNA.fa: $(GENCODE_DIR)/gencode.v32.transcripts.fa.gz selected_lncRNA
	zcat $< | fasta2oneline | tr "|" "\t" | filter_1col 6 <(cut -f 1 $^2) | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,10 | tab2fasta | fold > $@

#######################################
#
#	tipical final targets
#


#a tipical final target is
#someting.fa.tpx.tts_genom_coords.bed


#######################################
#
#	generic rules
#

%.header_added: %
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@
%.header_added.gz: %.gz
	(bawk -M $< | cut -f 2 | transpose; zcat $< ) | gzip > $@
%.xlsx: %.gz
	zcat $< | tab2xlsx > $@
%.xlsx: %
	cat $< | tab2xlsx > $@

%.slop_tpx.bed: %.bed
	bedtools sort -i $< | bedtools slop -l $(TPX_DISTANCE) -r $(TPX_DISTANCE) -i stdin -g $(GENCODE_DIR)/chrom.info > $@

%.merged.bed: %.bed
	bedtools sort < $< | bedtools merge -c 4 -o distinct > $@

%.mm9.bed: %.mm10.bed ../../local/share/data/mm10ToMm9.over.chain.gz
	liftOver $< $^2 $@ $@.unmap

##########################################3
#
#	region definition
#

cCRE.bed: /sto1/ref/bioinfotree/task/encode-screen/dataset/20201125/mm9-ES_E14.bed
	bawk '$$10=="PLS" || $$10=="dELS" || $$10=="pELS"' $< > $@

.bed: H3K27ac_selected.pre
	bawk '{print $$1~4}' $< > $@

%.bed.fa: $(GENCODE_DIR)/GRCh38.primary_assembly.genome.clean_id.fa %.bed
	bedtools getfasta -name -fi $< -bed $^2 -fo $@

##########################################3
#
#	genes association
#

#tss.enhancers.list: $(GENCODE_DIR)/primary_assembly.annotation.tss ../../local/share/data/1-s2.0-S0092867417311376-mmc3.csv
#	(bawk '$$1 {print $$1~4,"TSS"}' $<; bawk '$$1 && NR>3 {print $$1~4,"enhancer"}' $^2) | bsort -k1,1V -k2,2n > $@
#tss.enhancers.list.mm9: tss.enhancers.list ../../local/share/data/mm10ToMm9.over.chain.gz
#	liftOver $< $^2 $@ $@.unmap
#tss.enhancers.list.mm9.slop: tss.enhancers.list.mm9
#	bedtools sort -i $< | bedtools slop -l $(REGREGION_DISTANCE) -r $(REGREGION_DISTANCE) -i stdin -g $(BIOINFO_REFERENCE_ROOT)/gencode/dataset/mmusculus/M1/chrom.info | bedtools sort -i stdin > $@

tss.slop: $(GENCODE_DIR)/primary_assembly.annotation.tss
	bedtools sort -i $< | bedtools slop -l $(REGREGION_DISTANCE) -r $(REGREGION_DISTANCE) -i stdin -g $(GENCODE_DIR)/chrom.info | bedtools sort -i stdin > $@

%.regulated_genes.bed: %.bed tss.slop
	bedtools intersect -a $< -b $^2 -loj > $@

%.narrowPeak.mm9: %.narrowPeak.mm10.gz ../../local/share/data/mm10ToMm9.over.chain.gz
	liftOver <(zcat $< | cut -f -3,9) $^2 $@ $@.unmap

##########################################3
#
#	RNAseq
#

GEP.count.gz: /sto1/epigen/Maldotti_lncSMAD7/lavoro_lincRNA_Smad7/downstream_analysis/lncSmad7/Results/lncSmad7_rescue_Smad7exogenous_v7-gencodeM20/raw_counts.csv
	tr ";" "\t" < $< | grep_columns -k 1 NT_WT | bawk '$$2 > 0 || $$3>0 '> $@
GEP.count.ltmm.gz: GEP.count.gz
	R --slave --vanilla -e 'library(edgeR);\\
x <- read.table("$<", header=T,check.names=FALSE,row.names=1);\\
y <- DGEList(counts=x);\\
y <- calcNormFactors(y,method="TMM");\\
write.table(y$$samples,"$@.factors", sep="\t", quote=F, col.names=NA, row.names=T);\\
write.table(cpm(y, normalized.lib.sizes=TRUE, log=TRUE), "$@.tmp", sep="\t", quote=F, col.names=NA, row.names = T)'
	(echo -n "Geneid"; cat $@.tmp) | gzip > $@
	rm $@.tmp
RNA.exp.gz: GEP.count.ltmm.gz
	bawk 'BEGIN{print "GeneID","ltmm"} NR>1 {print $$1,($$3+$$2)/2}' $< | gzip > $@

.META: bulk_RNAseq.lncRNA_upInOld.H1hESC_localization.gz
	1	contrast
	2	GeneID
	3	logFC
	4	Pvalue
	5	Pvalue_adj
	6	significance
	7	biotype
	8	cytosol
	9	nucleus
	10	ratio2


###########################################
#
#	Acetilation
#

H3K27ac.encode_idr_conservative.narrowPeak.mm10.gz:
	wget -O $@ https://www.encodeproject.org/files/ENCFF512HTY/@@download/ENCFF512HTY.bed.gz
H3K27ac.ourV1.narrowPeak.gz: /sto1/epigen/Maldotti_lncSMAD7/lavoro_lincRNA_Smad7/ChIPseq/H3K27ac/peaks_narrow/dataset/v1/H3K27ac_WT_q0.05_fe2_peaks.narrowPeak
	cat $< | gzip > $@

H3K27ac.encode_idr_conservative.mm10.bw:
	wget -c -O $@ https://www.encodeproject.org/files/ENCFF583WVZ/@@download/ENCFF583WVZ.bigWig

%.H3K27ac.ourV1: %.bed ../../local/share/data/Qnorm500/H3K27ac_WT.500bp_quantile_normalized.bw
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) quay.io/biocontainers/ucsc-bigwigaverageoverbed:357--1 bash -c "cd $(PWD); bigWigAverageOverBed $^2 <(cut -f -4 $<) /dev/stdout -minMax" > $@
%.H3K27ac.encode_idr_conservative: %.bed H3K27ac.encode_idr_conservative.bw
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) quay.io/biocontainers/ucsc-bigwigaverageoverbed:357--1 bash -c "cd $(PWD); bigWigAverageOverBed $^2 <(cut -f -4 $<) /dev/stdout -minMax" > $@

.META: %.H3K27ac.ourV1 %.H3K27ac.encode_idr_conservative
	1	name	- name field from bed, which should be unique
	2	size 	- size of bed (sum of exon sizes)
	3	covered	- # bases within exons covered by bigWig
	4	sum	- sum of values over all bases covered
	5	mean0	- average over bases with non-covered bases counting as zeroes
	6	mean
	7	min
	8	max

##########################################3
#
#	triplexator
#


#%.tfo: %.fa
#	ssh epigen \
	sudo su - edoardo \
	/home/edoardo/src/triplexator/bin/triplexator $(TRIPLEXATOR_PARAM) -fm 1 -od . -o $@ -po -rm 2 -p 24 -ss $<
#%.tfo.log: %.tfo
#	@echo done
#%.tfo.summary: %.tfo
#	@echo done

%fa.tts: %fa
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) triplexator:v0.02 bash -c "cd $(PWD); triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 0     -o $@ -rm 2 -p $(CORES) -ds $<" 

%fa.tfo: %fa
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) triplexator:v0.02 bash -c "cd $(PWD); triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 0     -o $@ -rm 2 -p $(CORES) -ss $<"

%fa.tpx: ssRNA.fa %fa
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) triplexator:v0.02 bash -c "cd $(PWD); triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 0     -o $@ -rm 2 -p $(CORES) -ss $< -ds $^2" 
%fa.tpx_aln: ssRNA.fa %fa
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) triplexator:v0.02 bash -c "cd $(PWD); triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 1 -po -o $@ -rm 2 -p $(CORES) -ss $< -ds $^2" 

.META: *fa.tpx
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate

%.fa.tpx.around_validated_tfo: %.fa.tpx
	bawk '{named_tfo="no_tpx"; \
		if($$Duplex_ID!="no_tpx"){ \
			named_tfo="no_name";\
			if ($$TFO_start>=83 && $$TFO_start<=89) {named_tfo="86"} \
			if($$TFO_start>=45 && $$TFO_start<=51)  {named_tfo="48"} \
		} \
		print $$0,named_tfo} \
	' $< > $@

.META: *fa.around_validated_tfo
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	tfo_name

#########################################################
#
#	from tpx to bed (genome coords)
#

%.fa.tpx.tts_genom_coords.pre: %.fa.tpx %
	unhead $< | translate -r -a -f 4 $^2 4 > $@
%.fa.tpx.around_validated_tfo.tts_genom_coords.pre: %.fa.tpx.around_validated_tfo %
	unhead $< | translate -r -a -f 4 $^2 4 > $@

.META: *.tpx.tts_genom_coords.pre
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	region_chr
	15	region_b
	16	region_e

.META: *.tpx.tts_genom_coords
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	region_chr
	15	region_b
	16	region_e
	17	TTS_start_genom
	18	TTS_end_genom

.META: Srsf3_region.bed.fa.tpx.around_validated_tfo.tts_genom_coords
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	tfo_name
	15	region_chr
	16	region_b
	17	region_e
	18	TTS_start_genom
	19	TTS_end_genom

.META: *.tpx.around_validated_tfo.tts_genom_coords.pre *.tpx.around_validated_tfo.tts_genom_coords
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	tfo_name
	15	region_chr
	16	region_b
	17	region_e
	18	TTS_start_genom
	19	TTS_end_genom

%.tpx.tts_genom_coords: %.tpx.tts_genom_coords.pre
	bawk '{TTS_start_genom=$$TTS_start+$$region_b; TTS_end_genom=$$TTS_end+$$region_b; print $$0,TTS_start_genom,TTS_end_genom}' $< > $@
%.tpx.around_validated_tfo.tts_genom_coords: %.tpx.around_validated_tfo.tts_genom_coords.pre
	bawk '{TTS_start_genom=$$TTS_start+$$region_b; TTS_end_genom=$$TTS_end+$$region_b; print $$0,TTS_start_genom,TTS_end_genom}' $< > $@

%.tpx.tts_genom_coords.bed: %.tpx.tts_genom_coords
	bawk '{print $$region_chr,$$TTS_start_genom,$$TTS_end_genom, $$sequence_id "-" $$Duplex_ID, $$Score,$$Strand,$$sequence_id,$$TFO_start,$$TFO_end,$$Error_rate,$$Errors,$$Motif,$$Orientation,$$Guanine_rate,$$Duplex_ID,$$region_b,$$region_e}' $< \
	| id2count -a 4  >$@
%.tpx.around_validated_tfo.tts_genom_coords.bed: %.tpx.around_validated_tfo.tts_genom_coords
	bawk '{print $$region_chr,$$TTS_start_genom,$$TTS_end_genom, $$sequence_id "-" $$tfo_name "-" $$Duplex_ID, $$Score,$$Strand,$$sequence_id,$$TFO_start,$$TFO_end,$$Error_rate,$$Errors,$$Motif,$$Orientation,$$Guanine_rate,$$tfo_name,$$Duplex_ID,$$region_b,$$region_e}' $< \
	| id2count -a 4 >$@

.META: *.tpx.around_validated_tfo.tts_genom_coords.bed
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	tfo_name
	16	region_name
	17	region_b
	18	region_e

.META: *.tpx.tts_genom_coords.bed
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	region_name
	16	region_b
	17	region_e


.META: dba_results.contrast_Old_vs_Young.fdr_0.05.gain.atac_gain.bed.fa.tpx.tts_genom_coords.cutoffs.tpx.peaks.perc
	1	GeneID
	2	tpx_cutoff_10
	3	peaks_cutoff_10
	4	perc_peaks_cutoff_10
	5	tpx_cutoff_15
	6	peaks_cutoff_15
	7	perc_peaks_cutoff_15
	8	tpx_cutoff_18
	9	peaks_cutoff_18
	10	perc_peaks_cutoff_18


##########################################################################
#
#	repeat in tts
#

tts.simple_repeat: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed /sto1/ref/bioinfotree/task/ucsc-tracks/dataset/mm9/repeat_simple.extended.bed.gz
	cut -f -4 $< | bedtools intersect -a - -b $^2 -loj | cut -f 4,8 > $@

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.repeat_simple: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed /sto1/ref/bioinfotree/task/ucsc-tracks/dataset/mm9/repeat_rmsk.simple.bed.gz
	bedtools intersect -a $< -b $^2 -loj > $@
cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.repeat_noSimple: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed /sto1/ref/bioinfotree/task/ucsc-tracks/dataset/mm9/repeat_rmsk.noSimple.bed.gz
	bedtools intersect -a $< -b $^2 -loj > $@

.META: *.tpx.tts_genom_coords.bed.repeat_noSimple *.tpx.tts_genom_coords.bed.repeat_simple
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	region_name
	16	region_b
	17	region_e
	18	repeat_chr
	19	repeat_b
	20	repeat_e
	21	repeat_name

tts_simple_repeat_enric: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.repeat_simple
	parallel 'bawk '\''$$Score>{} {print $$1~3,$$repeat_name}'\'' $< | bedtools sort | bedtools merge -c 4 -o distinct | cut -f 4 | symbol_count -n | bsort -k2,2nr | bawk '\''$$2>=0.001'\'' | append_each_row -B {}' ::: 13 16 19 21 | tab2matrix -r repeat -t > $@



####################################################################
#
#	integration
#

cCRE.H3K27ac.%.peak_qscore: cCRE.bed H3K27ac.%.narrowPeak.gz
	bedtools intersect -a $< -b <(zcat $^2 | cut -f -3,9) -loj | cut -f 4,15 | bawk '$$2!="."' | bsort | stat_base -g -b > $@
cCRE.H3K27ac.%.peak_qscore: cCRE.bed H3K27ac.%.narrowPeak.mm9
	bedtools intersect -a $< -b $^2 -loj | cut -f 4,15 | bawk '$$2!="."' | bsort | stat_base -g -b > $@

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed cCRE.bed
	bedtools closest -a <(bedtools sort < $<) -b <(bedtools sort < $^2 | cut -f 1-4,10) -d | cut -f -14,19- > $@

.META: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	cCRE_b
	16	cCRE_e
	17	cCRE_name
	18	cCRE_type
	19	cCRE_TTS_distance


H3K27ac_selected.regulated_genes.RNAseq.gz: H3K27ac_selected.regulated_genes.bed RNAseq.logFC
	 cut -f 4,8,9 $< | translate -a -r -v -e NA $^2 2 | bsort | uniq | gzip > $@
	
.META: H3K27ac_selected.regulated_genes.RNAseq.gz
	1	peak_name	r1
	2	GeneID		Prdm14
	3	regreg_type	enhancer
	4	RNA_Rlnc_KO_logFC
	5	RNA_Rsmad_KO_logFC
	6	RNA_WT_KO_logFC

H3K27ac_selected.regulated_genes.RNAseq.with_tpx.gz: H3K27ac_selected.regulated_genes.RNAseq.header_added.gz H3K27ac_selected.bed.fa.tpx.around_validated_tfo.header_added
	zcat $< | translate -a -v -e no_tpx -d -j -f 4 <(sed s/Duplex_ID/peak_name/ $^2) 4 | gzip > $@

.META: H3K27ac_selected.regulated_genes.RNAseq.with_tpx.around_validated_tfo.gz
	1	regreg_chr	chr1
	2	regreg_b	13065030
	3	regreg_e	13075031
	4	gene		Prdm14
	5	regreg_type	enhancer
	6	peak_chr	chr1
	7	peak_b		13074607
	8	peak_e		13074757
	9	peak_name	r1
	10	Duplex_ID
	11	TFO_start
	12	TFO_end
	13	TTS_start
	14	TTS_end
	15	Score
	16	Error_rate
	17	Errors
	18	Motif
	19	Strand
	20	Orientation
	21	Guanine_rate
	22	K27ac_WT_KO	0.52299
	23	K27ac_R_KO	0.52889
	24	RNA_KO_WT	1.07522090383692
	25	RNA_KOSmad_WT	0.503384768453576
	26	KOTlncSmad_WT	0.570044409040709
	27	tfo_name	no_name|48|86

RNA_all: KO-WT.DGE_complete.txt KOTlncSmad-WT.DGE_complete.txt KOTSmad-WT.DGE_complete.txt
	cat $< | translate -r -a -v -e no_lncSmad7 $^2 1 | translate -r -a -v -e no_Smad7 $^3 1 | unhead > $@

.META: RNAselection RNA_all
	1	gene			Dnmt3a
	2	KO_WT_logFC 		-2.62382261477289
	3	KO_WT_CPM		7.16289031577841
	4	KO_WT_LR		1763.27176855942
	5	KO_WT_Pvalue		0
	6	KO_WT_FDR		0
	7	KOTlncSmad_WT_logFC	-2.62382261477289
	8	KOTlncSmad_WT_CPM	7.16289031577841
	9	KOTlncSmad_WT_LR	1763.27176855942
	10	KOTlncSmad_WT_Pvalue	0
	11	KOTlncSmad_WT_FDR	0
	12	KOTSmad_WT_logFC	-2.62382261477289
	13	KOTSmad_WT_CPM		7.16289031577841
	14	KOTSmad_WT_LR		1763.27176855942
	15	KOTSmad_WT_Pvalue	0
	16	KOTSmad_WT_FDR		0


genes_rnaseq_selected: /sto1/epigen/Maldotti_lncSMAD7/lavoro_lincRNA_Smad7/downstream_analysis/lncSmad7/Results/lncSmad7_rescue_Smad7exogenous_v5/target_lncSmad7_specific_ordered_WT_KO.txt
	cp $< $@

H3K27ac_selected.regulated_genes.RNAseq.with_tpx.around_validated_tfo.rnaseq_selected.gz: H3K27ac_selected.regulated_genes.RNAseq.with_tpx.around_validated_tfo.gz genes_rnaseq_selected
	zcat $< | translate -r -a -v -e no_rnaseq_selected <(bawk '{print $$1,"rnaseq_selected"}' $^2) 4 | gzip > $@

%.uniq.bed: %.bed
	bedtools sort < $< | uniq > $@


############################################
#
#	integration new
#


cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.ourV1.peak_qscore: cCRE.H3K27ac.ourV1.peak_qscore cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE
	translate -a -v -e NA -r $< 17 < $^2 > $@
cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac: cCRE.H3K27ac.encode_idr_conservative.peak_qscore cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.ourV1.peak_qscore
	translate -a -v -e NA -r $< 17 < $^2 > $@

.META: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	cCRE_b
	16	cCRE_e
	17	cCRE_name
	18	cCRE_type
	19	cCRE_TTS_distance
	20	H3K27ac_ourV1_qscore
	21	H3K27ac_encode_idr_conservative_qscore

cCRE.regulated_genes: cCRE.regulated_genes.bed
	bawk '{g=($$13+$$14)/2; if(d<0){d=c-g}; print $$4,$$15,$$16,int(g)}' $< > $@
cCRE.regulated_genes.RNA: RNA.exp.gz cCRE.regulated_genes
	translate -a -r -v -e NA <(zcat $<) 2 < $^2 > $@

.META: cCRE.regulated_genes.RNA
	1	cCRE	EM10E0431291
	2	GeneID	Gm37144
	3	gene_reference_poinit_type	TSS
	4	gene_reference_poinit	6285.5
	5	RNA_WT_ltmm	0.61955780234172797

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA: cCRE.regulated_genes.RNA cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac
	translate -a -d -j $< 17 < $^2 > $@

.META: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	cCRE_b
	16	cCRE_e
	17	cCRE_name
	18	GeneID	Gm37144
	19	gene_reference_poinit_type	TSS
	20	gene_reference_poinit	6285.5
	21	RNA_WT_ltmm	0.61955780234172797
	22	cCRE_type
	23	cCRE_TTS_distance
	24	H3K27ac_ourV1_qscore
	25	H3K27ac_encode_idr_conservative_qscore

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.pre: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA
	bawk '{t=($$TTS_start+$$TTS_end)/2; d=$$gene_reference_poinit-t; if(d<0){d=t-$$gene_reference_poinit}; print $$0,d}' $< > $@

.META: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.pre
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	cCRE_b
	16	cCRE_e
	17	cCRE_name
	18	GeneID	Gm37144
	19	gene_reference_poinit_type	TSS
	20	gene_reference_poinit	6285.5
	21	RNA_WT_ltmm	0.61955780234172797
	22	cCRE_type
	23	cCRE_TTS_distance
	24	H3K27ac_ourV1_qscore
	25	H3K27ac_encode_idr_conservative_qscore
	26	dist

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.pre
	bawk '{print $$TTS_name","$$GeneID","$$gene_reference_poinit_type,$$dist,$$0}' $< | find_best 1 2 -r | cut -f 3- > $@

.META: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.pre
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	cCRE_b
	16	cCRE_e
	17	cCRE_name
	18	GeneID	Gm37144
	19	gene_reference_poinit_type	TSS
	20	gene_reference_poinit	6285.5
	21	RNA_WT_ltmm	0.61955780234172797
	22	cCRE_type
	23	cCRE_TTS_distance
	24	H3K27ac_ourV1_qscore
	25	H3K27ac_encode_idr_conservative_qscore
	26	dist

.META: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	cCRE_b
	16	cCRE_e
	17	cCRE_name
	18	GeneID	Gm37144
	19	gene_reference_poinit_type	TSS
	20	gene_reference_poinit	6285.5
	21	RNA_WT_ltmm	0.61955780234172797
	22	cCRE_type
	23	cCRE_TTS_distance
	24	H3K27ac_ourV1_qscore
	25	H3K27ac_encode_idr_conservative_qscore
	26	dist
	27	repeatNoSimple
	28	repeatSimple

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.repeat_noSimple cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.repeat_simple
	cat $< \
	| translate -a -d -r <(bawk '{print $$TTS_name,$$repeat_name}' $^2 | tr ";" ",") 4 \
	| translate -a -d -r <(bawk '{print $$TTS_name,$$repeat_name}' $^3 | tr ";" ",") 4  > $@

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.fa cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat
	translate -a -r <(fasta2oneline $<) 4 < $^2 > $@

.META: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter.Malat1 cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter.Malat1.deoverlap cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.Malat1
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	cCRE_b
	16	cCRE_e
	17	cCRE_name
	18	GeneID	Gm37144
	19	gene_reference_poinit_type	TSS
	20	gene_reference_poinit	6285.5
	21	RNA_WT_ltmm	0.61955780234172797
	22	cCRE_type
	23	cCRE_TTS_distance
	24	H3K27ac_ourV1_qscore
	25	H3K27ac_encode_idr_conservative_qscore
	26	dist
	27	repeatNoSimple
	28	repeatSimple
	29	seq

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq
	bawk '$$RNA_WT_ltmm>1 &&  $$H3K27ac_ourV1_qscore!="NA" &&  $$H3K27ac_encode_idr_conservative_qscore!="NA"' $< > $@


cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.Malat1: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq
	bawk '$$sequence_id=="Malat1"' $< > $@
cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter.Malat1: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter
	bawk '$$sequence_id=="Malat1"' $< > $@

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter.Malat1.deoverlap.pre: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter.Malat1
	bawk '$$Score>=14 {print $$region_chr,$$TTS_start,$$TTS_end,$$TTS_name"/"$$Score,$$TFO_start,$$TFO_end}' $< | bedtools sort | bedtools merge -c 4 -o distinct -delim ";" | bawk '{print $$1";"$$2";"$$3,$$4}' | expandsets 2 | tr "/" "\t" | find_best 1 3 > $@
cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter.Malat1.deoverlap: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter.Malat1.deoverlap.pre cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter.Malat1
	filter_1col 4 <(cut -f 2 $<) < $^2 > $@

Malat1.deoverlap.tfo_profile: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.seq.filter.Malat1.deoverlap
	bawk '{print $$TFO_start,$$TFO_end,$$Score}' $< | perl -lane 'BEGIN{$$,="\t"} for($$F[0]..$$F[1]){print $$_,$$F[2]}' | symbol_count | bsort -n > $@

.META: Malat1.deoverlap.tfo_profile
	1	pos_in_transcript
	2	score
	3	tfo_count
	
Malat1.deoverlap.tfo_profile.png: Malat1.deoverlap.tfo_profile.header_added
	R --slave --vanilla -e 'library(ggplot2);\\
z<-read.delim("$<",header = T);\\
p<-ggplot(data=z,aes(x=pos_in_transcript,y=tfo_count,fill=score))+geom_bar(position="stack", stat="identity")+facet_wrap(~score, ncol=1, scales = "free_y");\\
ggsave("$@", width = 19, height = 10, dpi=100)'

%.fa.tpx_aln.oneline: %.fa.tpx_aln
	unhead $< | perl -pe '$$_=">$$_" if ($$.-1)%6==0' | perl -ne 'chomp; print "$$_\t"' | tr ">" "\n" | unhead > $@

.META: *fa.tpx_aln.oneline
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	alignment1
	15	alignment2
	16	alignment3
	17	alignment4

%.fa.tpx.tts_genom_coords.aln: %.fa.tpx.tts_genom_coords %.fa.tpx_aln.oneline
	bawk '{print 	$$Duplex_ID";"$$TTS_start";"$$TTS_end";"$$sequence_id";"$$TFO_start";"$$TFO_end, $$0}' $< \
	| translate -a -r <(bawk '{print \
			$$Duplex_ID";"$$TTS_start";"$$TTS_end";"$$sequence_id";"$$TFO_start";"$$TFO_end, $$Score, $$alignment1,$$alignment2,$$alignment3,$$alignment4}' $^2\
		| find_best 1 2 | cut -f 1,3- ) 1 \
	| cut -f 2- > $@

.META: *.fa.tpx.tts_genom_coords.aln
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	region_chr
	15	region_b
	16	region_e
	17	TTS_start_genom
	18	TTS_end_genom
	19	alignment1
	20	alignment2
	21	alignment3
	22	alignment4

%.EP300_encode.bed: %.bed ../igv_traks/EP300_BruceES.idr_condservative.mm9.named.narrowpeak.gz
	bedtools intersect -loj -a <(cut -f -4 $< | bedtools sort | uniq) -b $^2 | bawk '{if($$13=="."){$$13=-1}; print}' | find_best -H 4 13 > $@

%.EP300_Song.bed: %.bed ../igv_traks/EP300_Song.mm9.bed.gz
	bedtools intersect -loj -a <(cut -f -4 $< | bedtools sort | uniq) -b $^2 | bawk '{if($$9=="."){$$9=-1}; print}' | find_best -H 4 9 > $@

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.EP300_encode.bed cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.EP300_Song.bed
	translate -a -r <(bawk '{print $$4,$$13}' $<) 4 < $^2 | translate -a -r <(bawk '{print $$4,$$9}' $^3) 4 > $@

.META: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300 
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	cCRE_b
	16	cCRE_e
	17	cCRE_name
	18	GeneID	Gm37144
	19	gene_reference_poinit_type	TSS
	20	gene_reference_poinit	6285.5
	21	RNA_WT_ltmm	0.61955780234172797
	22	cCRE_type
	23	cCRE_TTS_distance
	24	H3K27ac_ourV1_qscore
	25	H3K27ac_encode_idr_conservative_qscore
	26	dist
	27	repeatNoSimple
	28	repeatSimple
	29	EP300_encode_qvalue
	30	EP300_Song_score

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.header_added cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.aln.header_added
	bawk '{print $$region_chr";"$$TTS_start";"$$TTS_end";"$$sequence_id";"$$TFO_start";"$$TFO_end, $$0}' $< | translate -a -r <(\
		bawk '{print $$region_chr ";" $$TTS_start_genom ";" $$TTS_end_genom";"$$sequence_id";"$$TFO_start";"$$TFO_end,$$Score,$$19~22}' $^2 \
		| find_best -H 1 2 | cut -f 1,3- | perl -pe 's/_genom//g if $$.==1') 1 \
	| cut -f 2- > $@
cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP: PARalyzer.RIP.header_added cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln
	translate -a -f 2 <(sed 's/GeneID/sequence_id/' < $<) 7 < $^2 > $@

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP.filt: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP
	bawk '$$RNA_WT_ltmm>1 && ($$EP300_Song_score!=-1 || $$EP300_encode_qvalue!=-1) && ($$H3K27ac_ourV1_qscore!="NA" || $$H3K27ac_encode_idr_conservative_qscore!="NA")' $< > $@

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP.filt.5930430L01Rik cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP.filt.B230208H11Rik cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP.filt.Malat1 cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP.filt.Mir155hg cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP.filt.1700100I10Rik: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP.filt.%: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.cCRE.H3K27ac.RNA.dist.repeat.P300.aln.PACRLIP.RIP.filt
	bawk 'NR==1 || $$sequence_id=="$*"' $< > $@

pre_tpx-closedYoung_openOld_diffBindOld.bed: /sto1/epigen/Polignano_P300senescence/complex-analysis/dataset/v3/pre_tpx-closedYoung_openOld_diffBindOld.bed
	link_install $< $@

CORES=10
RGTDATA_DIR=/sto1/ref/bioinfotree/task/rgtdata/0.13.1
DOCKER_GROUP=1000


%/index.html: pre_tpx-closedYoung_openOld_diffBindOld.bed %
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $$PWD:$$PWD -v $(RGTDATA_DIR):/opt/rgtdata -e RGTDATA=/opt/rgtdata rgt:0.13.1 bash -c "cd $$PWD; rgt-TDF regiontest -r $(word 2,$^) -bed $(word 1,$^) -rn $(word 2,$^) -organism hg38 -o $(word 2,$^)TDF_1L_fullnew_l10e20E1 -ccf 1 -n 1000 -mp $(CORES) -l 10 -e 20 -par L_-1_-E_-1"

%/default_index.html: pre_tpx-closedYoung_openOld_diffBindOld.bed %
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $$PWD:$$PWD -v $(RGTDATA_DIR):/opt/rgtdata -e RGTDATA=/opt/rgtdata rgt:0.13.1 bash -c "cd $$PWD; rgt-TDF regiontest -r $(word 2,$^) -bed $(word 1,$^) -rn $(word 2,$^) -organism hg38 -o $(word 2,$^)TDF_default -n 1000 -mp $(CORES)"


TDF_ChIRP_01L_full_l10e20E1: pre_tpx-closedYoung_openOld_diffBindOld.bed ssRNA_filt1.fa
	docker run -u `id -u`:$(DOCKER_GORUP) --rm -v $$PWD:$$PWD -v $(RGTDATA_DIR):/opt/rgtdata -e RGTDATA=/opt/rgtdata rgt:0.13.1 bash -c "cd $$PWD; rgt-TDF regiontest -r $(word 2,$^) -bed $(word 1,$^) -rn ssRNA_filt1 -organism hg38 -o TDF_01L_l10e20E1 -ccf 0.1 -n 1000 -mp $(CORES) -l 10 -e 20 -par L_-1_-E_-1"


%.gat.done: /sto1/epigen/Polignano_P300senescence/local/share/rules/rules-gat.mk
	mkdir $*.faTDF_1L_fullnew_l10e20E1/$*.fa/gat; \
	cd $*.faTDF_1L_fullnew_l10e20E1/$*.fa/gat; \
	ln -s $< rules.mk; touch makefile; \
	remake NCPU=6 $*.gat_runs_done; \
	cd /sto1/epigen/Polignano_P300senescence/tpx_analysis/dataset/TDF; \
	touch $@
ssRNAs_selected.DBD_target_regions.best_score.bed: AC073283.2.faTDF_1L_fullnew_l10e20E1/AC073283.2.fa/grafo/AC073283.2.fa_DBD_target_regions.best_score.bed LINC00973.faTDF_1L_fullnew_l10e20E1/LINC00973.fa/grafo/LINC00973.fa_DBD_target_regions.best_score.bed LINC02154.faTDF_1L_fullnew_l10e20E1/LINC02154.fa/grafo/LINC02154.fa_DBD_target_regions.best_score.bed LINC02575.faTDF_1L_fullnew_l10e20E1/LINC02575.fa/grafo/LINC02575.fa_DBD_target_regions.best_score.bed PURPL.faTDF_1L_fullnew_l10e20E1/PURPL.fa/grafo/PURPL.fa_DBD_target_regions.best_score.bed AC046195.1.faTDF_1L_fullnew_l10e20E1/AC046195.1.fa/grafo/AC046195.1.fa_DBD_target_regions.best_score.bed INKA2-AS1.faTDF_1L_fullnew_l10e20E1/INKA2-AS1.fa/grafo/INKA2-AS1.fa_DBD_target_regions.best_score.bed LINC01638.faTDF_1L_fullnew_l10e20E1/LINC01638.fa/grafo/LINC01638.fa_DBD_target_regions.best_score.bed MEG3.faTDF_1L_fullnew_l10e20E1/MEG3.fa/grafo/MEG3.fa_DBD_target_regions.best_score.bed MYOSLID.faTDF_1L_fullnew_l10e20E1/MYOSLID.fa/grafo/MYOSLID.fa_DBD_target_regions.best_score.bed SFTA1P.faTDF_1L_fullnew_l10e20E1/SFTA1P.fa/grafo/SFTA1P.fa_DBD_target_regions.best_score.bed PCAT6.faTDF_1L_fullnew_l10e20E1/PCAT6.fa/grafo/PCAT6.fa_DBD_target_regions.best_score.bed
	cat $< $^2 $^3 $^4 $^5 $^6 $^7 $^8 $^9 $^10 $^11 $^12 | sed 's/.fa$$//' > $@
ssRNAs_selected.DBD_target_regions.best_score.cut10.grafo: ssRNAs_selected.DBD_target_regions.best_score.bed
	bawk '($$5>=10) {print $$1":"$$2"-"$$3,$$4~5,$$7}' $< > $@
ssRNAs_selected.DBD_target_regions.best_score.grafo: ssRNAs_selected.DBD_target_regions.best_score.bed
	bawk '{print $$1":"$$2"-"$$3,$$4~5,$$7}' $< > $@

motifs/%: %.bed
	conda activate bit_chipseq_0.1;\
        findMotifsGenome.pl $< hg38 $@ -size given -p $(CORES)
gat.repClass.matrix.transposed: AC073283.2.faTDF_1L_fullnew_l10e20E1/AC073283.2.fa/gat/gat.repClass.matrix.transposed LINC00973.faTDF_1L_fullnew_l10e20E1/LINC00973.fa/gat/gat.repClass.matrix.transposed LINC02154.faTDF_1L_fullnew_l10e20E1/LINC02154.fa/gat/gat.repClass.matrix.transposed LINC02575.faTDF_1L_fullnew_l10e20E1/LINC02575.fa/gat/gat.repClass.matrix.transposed PURPL.faTDF_1L_fullnew_l10e20E1/PURPL.fa/gat/gat.repClass.matrix.transposed AC046195.1.faTDF_1L_fullnew_l10e20E1/AC046195.1.fa/gat/gat.repClass.matrix.transposed INKA2-AS1.faTDF_1L_fullnew_l10e20E1/INKA2-AS1.fa/gat/gat.repClass.matrix.transposed LINC01638.faTDF_1L_fullnew_l10e20E1/LINC01638.fa/gat/gat.repClass.matrix.transposed MEG3.faTDF_1L_fullnew_l10e20E1/MEG3.fa/gat/gat.repClass.matrix.transposed MYOSLID.faTDF_1L_fullnew_l10e20E1/MYOSLID.fa/gat/gat.repClass.matrix.transposed SFTA1P.faTDF_1L_fullnew_l10e20E1/SFTA1P.fa/gat/gat.repClass.matrix.transposed PCAT6.faTDF_1L_fullnew_l10e20E1/PCAT6.fa/gat/gat.repClass.matrix.transposed
	cat $< $^2 $^3 $^4 $^5 $^6 $^7 $^8 $^9 $^10 $^11 $^12 > $@
gat.repFamily.matrix.transposed: AC073283.2.faTDF_1L_fullnew_l10e20E1/AC073283.2.fa/gat/gat.repFamily.matrix.transposed LINC00973.faTDF_1L_fullnew_l10e20E1/LINC00973.fa/gat/gat.repFamily.matrix.transposed LINC02154.faTDF_1L_fullnew_l10e20E1/LINC02154.fa/gat/gat.repFamily.matrix.transposed LINC02575.faTDF_1L_fullnew_l10e20E1/LINC02575.fa/gat/gat.repFamily.matrix.transposed PURPL.faTDF_1L_fullnew_l10e20E1/PURPL.fa/gat/gat.repFamily.matrix.transposed AC046195.1.faTDF_1L_fullnew_l10e20E1/AC046195.1.fa/gat/gat.repFamily.matrix.transposed INKA2-AS1.faTDF_1L_fullnew_l10e20E1/INKA2-AS1.fa/gat/gat.repFamily.matrix.transposed LINC01638.faTDF_1L_fullnew_l10e20E1/LINC01638.fa/gat/gat.repFamily.matrix.transposed MEG3.faTDF_1L_fullnew_l10e20E1/MEG3.fa/gat/gat.repFamily.matrix.transposed MYOSLID.faTDF_1L_fullnew_l10e20E1/MYOSLID.fa/gat/gat.repFamily.matrix.transposed SFTA1P.faTDF_1L_fullnew_l10e20E1/SFTA1P.fa/gat/gat.repFamily.matrix.transposed PCAT6.faTDF_1L_fullnew_l10e20E1/PCAT6.fa/gat/gat.repFamily.matrix.transposed
	cat $< $^2 $^3 $^4 $^5 $^6 $^7 $^8 $^9 $^10 $^11 $^12 > $^13
motifs/%.target_regions.cut12.cln/knownResults.qVal05.txt: motifs/%.target_regions.cut12.cln/knownResults.txt
	bawk '$$5<=0.05' $< | awk -F "(" '{print "$*""\t"$$1"\t"$$0}' > $@
motifs.knownResults.qVal05.txt: motifs/AC073283.2.target_regions.cut12.cln/knownResults.qVal05.txt motifs/LINC02154.target_regions.cut12.cln/knownResults.qVal05.txt motifs/LINC02575.target_regions.cut12.cln/knownResults.qVal05.txt motifs/PURPL.target_regions.cut12.cln/knownResults.qVal05.txt motifs/AC046195.1.target_regions.cut12.cln/knownResults.qVal05.txt motifs/INKA2-AS1.target_regions.cut12.cln/knownResults.qVal05.txt motifs/SFTA1P.target_regions.cut12.cln/knownResults.qVal05.txt
	cat $< motifs/LINC00973.target_regions.cut12.cln/knownResults.qVal05.txt $^2 $^3 $^4 $^5 $^6 motifs/LINC01638.target_regions.cut12.cln/knownResults.qVal05.txt motifs/MEG3.target_regions.cut12.cln/knownResults.qVal05.txt motifs/MYOSLID.target_regions.cut12.cln/knownResults.qVal05.txt $^7 motifs/PCAT6.target_regions.cut12.cln/knownResults.qVal05.txt > $^8
motifs.knownResults.qVal05.txt: motifs/AC073283.2.target_regions.cut12.cln/knownResults.qVal05.txt motifs/LINC02154.target_regions.cut12.cln/knownResults.qVal05.txt motifs/LINC02575.target_regions.cut12.cln/knownResults.qVal05.txt motifs/PURPL.target_regions.cut12.cln/knownResults.qVal05.txt motifs/AC046195.1.target_regions.cut12.cln/knownResults.qVal05.txt motifs/INKA2-AS1.target_regions.cut12.cln/knownResults.qVal05.txt motifs/SFTA1P.target_regions.cut12.cln/knownResults.qVal05.txt
	cat $< motifs/LINC00973.target_regions.cut12.cln/knownResults.qVal05.txt $^2 $^3 $^4 $^5 $^6 motifs/LINC01638.target_regions.cut12.cln/knownResults.qVal05.txt motifs/MEG3.target_regions.cut12.cln/knownResults.qVal05.txt motifs/MYOSLID.target_regions.cut12.cln/knownResults.qVal05.txt $^7 motifs/PCAT6.target_regions.cut12.cln/knownResults.qVal05.txt > $^8
motifs.target_regions.cut12.cln: motifs/AC073283.2.target_regions.cut12.cln/knownResults.qVal05.txt motifs/LINC02154.target_regions.cut12.cln/knownResults.qVal05.txt motifs/LINC02575.target_regions.cut12.cln/knownResults.qVal05.txt motifs/PURPL.target_regions.cut12.cln/knownResults.qVal05.txt motifs/AC046195.1.target_regions.cut12.cln/knownResults.qVal05.txt motifs/INKA2-AS1.target_regions.cut12.cln/knownResults.qVal05.txt motifs/LINC01638.target_regions.cut12.cln/knownResults.qVal05.txt motifs/MEG3.target_regions.cut12.cln/knownResults.qVal05.txt motifs/MYOSLID.target_regions.cut12.cln/knownResults.qVal05.txt motifs/SFTA1P.target_regions.cut12.cln/knownResults.qVal05.txt motifs/PCAT6.target_regions.cut12.cln/knownResults.qVal05.txt
	cat $< $^2 $^3 $^4 $^5 $^6 $^7 $^8 $^9 $^10 $^11 > $^12

merged_similarLnc_target_regions.cut12.cln.bed: AC046195.1.target_regions.cut12.cln.bed PCAT6.target_regions.cut12.cln.bed SFTA1P.target_regions.cut12.cln.bed MEG3.target_regions.cut12.cln.bed AC073283.2.target_regions.cut12.cln.bed
	cat $< $^2 $^3 $^4 $^5 | bedtools sort -i stdin | bedtools merge -i stdin > $@

gat.rep%.matrix_done: ssRNAs_selected
	while read i; do cd "$$i".faTDF_1L_fullnew_l10e20E1/"$$i".fa/gat; remake gat.rep$*.matrix; cd -; done < $<

gat.rep%.header: gat.rep%.matrix_done
	cat *.faTDF_1L_fullnew_l10e20E1/*.fa/gat/gat.rep$*.matrix | cut -f1 | sort | uniq > $@

gat.rep%.matrix: gat.rep%.header AC073283.2.faTDF_1L_fullnew_l10e20E1/AC073283.2.fa/gat/gat.rep%.matrix LINC00973.faTDF_1L_fullnew_l10e20E1/LINC00973.fa/gat/gat.rep%.matrix LINC02154.faTDF_1L_fullnew_l10e20E1/LINC02154.fa/gat/gat.rep%.matrix LINC02575.faTDF_1L_fullnew_l10e20E1/LINC02575.fa/gat/gat.rep%.matrix PURPL.faTDF_1L_fullnew_l10e20E1/PURPL.fa/gat/gat.rep%.matrix AC046195.1.faTDF_1L_fullnew_l10e20E1/AC046195.1.fa/gat/gat.rep%.matrix INKA2-AS1.faTDF_1L_fullnew_l10e20E1/INKA2-AS1.fa/gat/gat.rep%.matrix LINC01638.faTDF_1L_fullnew_l10e20E1/LINC01638.fa/gat/gat.rep%.matrix MEG3.faTDF_1L_fullnew_l10e20E1/MEG3.fa/gat/gat.rep%.matrix MYOSLID.faTDF_1L_fullnew_l10e20E1/MYOSLID.fa/gat/gat.rep%.matrix SFTA1P.faTDF_1L_fullnew_l10e20E1/SFTA1P.fa/gat/gat.rep%.matrix PCAT6.faTDF_1L_fullnew_l10e20E1/PCAT6.fa/gat/gat.rep%.matrix
	cat $< | translate -a -r -v -e 0 <(cat $^2) 1 | translate -a -r -v -e 0 <(cat $^3 ) 1 | translate -a -r -v -e 0 <(cat $^4) 1 | translate -a -r -v -e 0 <(cat $^5) 1 | translate -a -r -v -e 0 <(cat $^6) 1 | translate -a -r -v -e 0 <(cat $^7) 1 | translate -a -r -v -e 0 <(cat $^8) 1 | translate -a -r -v -e 0 <(cat $^9) 1 | translate -a -r -v -e 0 <(cat $^10) 1 | translate -a -r -v -e 0 <(cat $^11) 1 | translate -a -r -v -e 0 <(cat $^12) 1 | translate -a -r -v -e 0 <(cat $^13) 1 | transpose > $@

gat.rep%.matrix.3: gat.rep%.matrix
	unhead $< | tr "_" "\t" > $@
ssRNAs_selected.DBD_target_regions.best_score.pre_NN.bed: ssRNAs_selected.DBD_target_regions.best_score.bed
	bawk '{print $$1~3 }' $< | sort | uniq > $^2
ssRNA.fa.tfo.2plot: ssRNA.fa.tfo
	bawk '{print $$1~3, $$2"-"$$3, $$4}' $< > $@
lncRNASmad7_full.fa.tfo.flt: lncRNASmad7_full.fa.tfo
	cat $< | bawk '{print $$1,$$2,$$3,$$2"-"$$3,$$4}' > $^2
