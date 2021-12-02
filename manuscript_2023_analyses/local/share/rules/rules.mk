##########################
#
#     General param
#

GENCODE_DIR=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/hsapiens/32
NCPU?=28


###########################
# 
#      Docker param
#

DOCKER_GROUP?=$(shell stat i-c %g $(PWD))
DOCKER_DATA_DIR?=$(PRJ_ROOT)



#######################################
#
#     Generic rules
#

%.header_added.bed: %.bed
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@
%.header_added: %
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@
%.header_added.gz: %.gz
	(bawk -M $< | cut -f 2 | transpose; zcat $< ) | gzip > $@
%.xlsx: %.gz
	zcat $< | tab2xlsx > $@
%.xlsx: %
	cat $< | tab2xlsx > $@
%.slop_tpx.bed: %.bed $(GENCODE_DIR)/chrom.info
	bedtools sort -i $< | bedtools slop -l $(TPX_DISTANCE) -r $(TPX_DISTANCE) -i stdin -g $^2 > $@
%.merged.bed: %.bed
	bedtools sort < $< | bedtools merge -c 4 -o distinct > $@
%.bb: %.bed $(GENCODE_DIR)/chrom.info
	bedToBigBed $< <(unhead $^2) $@
%.complement.bed: %.bed chrom.info.sorted
	bedtools complement -i <(bsort -k1,1V -k2,2n $<) -g $^2 > $@
chrom.info.sorted: $(GENCODE_DIR)/chrom.info
	unhead $< | bsort -k1,1V | grep -v _ > $@


######################################
#
#    LongTarget
#

cCRE.bed: $(BIOINFO_REFERENCE_ROOT)/encode-screen/dataset/v13/hg38-ccREs.bed
	ln -s $< $@

%.bed.fa: $(GENCODE_DIR)/GRCh38.primary_assembly.genome.fa %.bed
	bedtools getfasta -name -fi $< -bed $^2 -fo $@

%.mask.fa: $(GENCODE_DIR)/GRCh38.primary_assembly.genome.fa %.bed
	bedtools maskfasta -fi $< -fo $@ -bed $^2

chrs_mask/chr%: cCRE.complement.mask.fa
	mkdir -p `dirname $@`
	split_fasta -e -p $< 1 `dirname $@`
	for i in `dirname $@`/chr*; do perl -lpe '$$_.="|NA" if $$.==1' < $$i > $$i.tmp; mv $$i.tmp $$i; done

#all_lncRNA.LongTarget_mask.done: selected_ssRNA_id
#	for ssRNA in selected_ssRNA_id; do mkdir -p $$ssRNA; cd $$ssRNA; echo `seq 1 22` X Y M | tr " " "\n" | parallel -j 22 'LongTarget -f1 chr{} -f2 ssRNA.fa -ni 10 -ds 8 -lg 10'

%/ssRNA.fa: longest_transcripts.fa
	mkdir -p `dirname $@`
	get_fasta -i $* < $^ > $@


CHR?=chr1
#%/$(CHR)--$(CHR).t-TFOsorted: chrs_mask/$(CHR) chrs_mask/ssRNA.fa
#        cd chrs_mask;\
#        LongTarget -f1 $(CHR) -f2 ssRNA.fa -ni 10 -ds 8 -lg 10

%.ALL_LongTarget_mask.done: %/ssRNA.fa chrs_mask/chr1
	cd $*;\
	echo `seq 1 22` X Y M | tr " " "\n" | parallel -j 22 'LongTarget -f1 ../chrs_mask/chr{} -f2 ssRNA.fa -ni 10 -ds 8 -lg 10'
	touch $@



##################################
#
#    Plot fisher
#

ssRNA=AC004862.1 AC005520.5 AC008429.1 AC009158.1 AC009264.1 AC009299.3 AC012668.3 AC013652.1 AC015908.2 AC016590.1 AC019294.2 AC021752.1 AC025171.1 AC073130.1 AC073529.1 AC087477.2 AC087521.2 AC092296.2 AC097376.3 AC099520.1 AC103702.2 AC104461.1 AC105383.1 AC105450.1 AC114316.1 AC138627.1 AC139099.1 AC244502.1 ADAMTS9-AS1 AF165147.1 AL021396.1 AL031056.1 AL138828.1 AL355990.1 AL365259.1 AL390198.1 AL591368.1 AP000704.1 AP005263.1 C1RL-AS1 C5orf64 CALML3-AS1 CASC19 CHRM3-AS2 CLCA4-AS1 DNAH17-AS1 DNAJC27-AS1 DNM3OS EMX2OS EPHA1-AS1 FAM106A FAM230J FAM41C FAM87A GORAB-AS1 GTSE1-DT IGFL2-AS1 LACTB2-AS1 LEF1-AS1 LINC00323 LINC00467 LINC00616 LINC00837 LINC00877 LINC01127 LINC01378 LINC01471 LINC01483 LINC01531 LINC01583 LINC01606 LINC01828 LINC01829 LINC01934 LINC02068 LINC02112 LINC02389 LINC02453 LINC02762 MAP3K20-AS1 MBNL1-AS1 MEG3 MIR22HG MIR3976HG MYCBP2-AS1 MZF1-AS1 NAGPA-AS1 NEAT1 OTX2-AS1 PRDM16-DT SAMMSON SATB2-AS1 SGMS1-AS1 SMC2-AS1 SMIM2-AS1 SPATA41 TBX5-AS1 TCL6 TEX26-AS1 TMEM132D-AS1 TMEM99 VLDLR-AS1
COND=HEP_H9 NPC_H9 MESEND_H1 MESENC_H1


tpx_analysis.fisher_select_cutoff.ALLconditions.gz:
	matrix_reduce -t 'tpx_analysis/*/cCRE.tpx.best.complete.*_neg_pos.fisher_select_cutoff' | tr ";" "\t" | cut -f 1,2,4- | gzip > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.best: tpx_analysis.fisher_select_cutoff.ALLconditions.gz
	zcat $< | find_best -r 1:2:3 10 > $@

.META: tpx_analysis.fisher_select_cutoff.ALLconditions.gz tpx_analysis.fisher_select_cutoff.ALLconditions.best tpx_analysis.fisher_select_cutoff.ALLconditions.best.selected_ssRNA_id
	1	lncRNA			AC004862.1
	2	condition		MESENC_H1
	3	cCRE			pELS
	4	tpx_score		0.0
	5	greater_positive	24223
	6	lower_positive		0
	7	greater_negative	102094
	8	lower_negative		0
	9	oddsratio		nan
	10	pvalue			1
	11	TPXcCRE_score

#tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix: tpx_analysis.fisher_select_cutoff.ALLconditions.best.selected_ssRNA_id
#	bawk '$$cCRE=="$*" {$$pvalue=$$pvalue+2.26261e-288; print $$lncRNA,$$condition,-log($$pvalue)/log(10)}' $< | tab2matrix -r GeneID > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.matrix: tpx_analysis.fisher_select_cutoff.ALLconditions.best
	bawk '{$$pvalue=$$pvalue+2.26261e-288; print $$lncRNA,$$condition";"$$cCRE,-log($$pvalue)/log(10)}' $< | tab2matrix -r GeneID > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix.heatmap_euclid.pdf: tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix
	heatmap --pdf -e 13 -w 3 -n -c 12 -r 4 $< $@
tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix.heatmap_pearson.pdf: tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix
	heatmap --pdf -e 13 -w 3 -n -c 12 -r 4 -d correlation $< $@

GEP.count.exp_filter.ltmm.metadata.exp_genes_condition.matrix.selected_lncRNA.heatmap.pdf: GEP.count.exp_filter.ltmm.metadata.exp_genes_condition.matrix.selected_lncRNA.header_added
	heatmap -c 12 -e 13 -w 3 --pdf $< $@

all_lncRNA.ALLconditions:
	ls precomputed_tpx/* | perl -lne 'BEGIN{@lines=("NPC\tH9","HEP\tH9","MESEND\tH1","MESENC\tH1")} s|.*/||; s|-cCRE.*||; $$r=$$_; for(@lines){print "$$r\t$$_"}' > $@

%/cCRE.tpx.best.complete.NPC_H9_neg_pos.assoc_genes.gz: %/cCRE.tpx.best.complete.NPC_H9_neg_pos.gz %/cCRE.tpx.tts_genom_coords.assoc_genes cCRE.bed
	zcat $< | translate -a -d -j -v -e NA <(bawk '$$5!="." && $$5 {print $$2,$$3,$$4,$$5}' $^2 | sort | uniq | translate <(cut -f 4,5 $^3) 1) 2 | gzip > $@

.META: cCRE.tpx.best.complete.NPC_H9_neg_pos.assoc_genes.gz
	1	ssRNA	AC007098.1
	2	cCRE	EH38E0065969
	3	genehancer_region	NA
	4	assoc_gene	NA
	5	assoc_entrezGene	0
	6	tpx_score
	7	cCRE_type	dELS
	8	pos_neg	neg

tpx_analysis/AC007098.1/cCRE.tpx.best.complete.NPC_H9_neg_pos.assoc_genes.header_added.universe: tpx_analysis/AC007098.1/cCRE.tpx.best.complete.NPC_H9_neg_pos.assoc_genes.header_added.gz
	bawk '$$pos_neg=="pos" && $$assoc_entrezGene!="NA" {print $$assoc_entrezGene}' $< > $@
tpx_analysis/AC007098.1/cCRE.tpx.best.complete.NPC_H9_neg_pos.assoc_genes.header_added.selected_11: tpx_analysis/AC007098.1/cCRE.tpx.best.complete.NPC_H9_neg_pos.assoc_genes.header_added.gz
	bawk '$$pos_neg=="pos" && $$assoc_entrezGene!="NA" && $$tpx_score>=11 {print $$assoc_entrezGene}' $< > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.best.selected_ssRNA_id: selected_ssRNA_id tpx_analysis.fisher_select_cutoff.ALLconditions.best
	filter_1col 2 $< < $^2 > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix: tpx_analysis.fisher_select_cutoff.ALLconditions.best.selected_ssRNA_id
	bawk '$$cCRE=="$*" {$$pvalue=$$pvalue+2.26261e-288; if($$oddsratio>=1) {print $$condition,$$lncRNA,-log($$pvalue)/log(10)}else{print $$condition,$$lncRNA,log($$pvalue)/log(10)}}' $< | tab2matrix -r GeneID > $@
tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix.oddsratio: tpx_analysis.fisher_select_cutoff.ALLconditions.best.selected_ssRNA_id
	bawk '$$cCRE=="$*" {print $$condition,$$lncRNA,$$oddsratio-1}' $< | tab2matrix -r GeneID > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.best%.length_dist: tpx_analysis.fisher_select_cutoff.ALLconditions.best%
	cut -f1,3,4 $< | symbol_count > $@

.META: tpx_analysis.fisher_select_cutoff.ALLconditions.best%.length_dist
	1	condition
	2	cCRE
	3	tpx_score
	4	count


tpx_analysis.fisher_select_cutoff.ALLconditions.selected_ssRNA_id: tpx_analysis.fisher_select_cutoff.ALLconditions.gz selected_ssRNA_id
	zcat $< | filter_1col 2 $^2 > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.pvalue_correct.gz: tpx_analysis.fisher_select_cutoff.ALLconditions.gz
	zcat $< | pvalue_correct -a -c 10 | gzip > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.best.selected_ssRNA_id.GEPcount: tpx_analysis.fisher_select_cutoff.ALLconditions.best.selected_ssRNA_id GEP.count.exp_filter.ltmm.metadata.exp_genes_condition.matrix.selected_lncRNA.header_added
	bawk '{print $$2";"$$1,$$3,$$11}' $< | translate -a <(cut -f1,3- $^2 | matrix2tab | bawk '{print $$1";"$$2,$$3}') 1 | tr ";" "\t" > $@

.META: tpx_analysis.fisher_select_cutoff.ALLconditions.best.selected_ssRNA_id.GEPcount
	1	GeneID
	2	condition	HEP_H9
	3	exp	2.67169721233
	4	cCRE	PLS
	5	TPXcCRE_score	-131.317235


tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix.exp_corr: GEP.count.exp_filter.ltmm.metadata.exp_genes_condition.matrix.selected_lncRNA.header_added tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix
	./corr_plot.R $< < $^2 | bawk '{print $$5,$$1~4}'  > $@



#####################################
#
#	LongTarget Stability
#
cCRE.tpx.best.complete.fisher_select_cutoff.ALLcond.best.selected_ssRNA_id: $(addprefix tpx_analysis/, $(addsuffix /cCRE.tpx.best.complete.fisher_select_cutoff.ALLcond.best,$(ssRNA)))
	matrix_reduce -t 'tpx_analysis/*/cCRE.tpx.best.complete.fisher_select_cutoff.ALLcond.best' > $@	

tpx_analysis/%/cCRE.tpx.best.complete.fisher_select_cutoff.ALLcond.best:
	matrix_reduce -t 'tpx_analysis/$*/cCRE.tpx.best.complete.*_neg_pos.fisher_select_cutoff' | cut -f1,3- | find_best -a 1:2 10 > $@

.META: cCRE.tpx.best.complete.fisher_select_cutoff.ALLcond.best.selected_ssRNA_id
	1	GeneID	AC016590.1
	2	condition	MESENC_H1
	3	cCRE	dELS
	4	score	10.0
	5	greater_pos	6981
	6	lower_pos	12839
	7	greater_neg	5843
	8	lower_neg	6037
	9	odds_ratio	0.561787
	10	pvalue	4.81687e-132
	11	TPXcCRE_score	-131.317235

#####################################

ALL.cCRE.tpx.best.complete.gz: cCRE.bed
	zcat tpx_analysis/*/cCRE.tpx.best.complete.gz | tab2matrix -t | translate -a -v -e pos_lncRNA <(cut -f 4,5 $<) 1 | gzip > $@

ALL.cCRE.tpx.best.complete.neg_pos.fisher_select_cutoff.best:
	matrix_reduce -t 'tpx_analysis/*/cCRE.tpx.best.complete.*.fisher_select_cutoff' | find_best -r 2 10 > $@

tpx_analysis/Z85996.1/cCRE.tpx.best.complete.Z85996.1-neg_pos.best_cut.gz: ALL.cCRE.tpx.best.complete.neg_pos.fisher_select_cutoff.best tpx_analysis/Z85996.1/cCRE.tpx.best.complete.Z85996.1-neg_pos.gz
	bawk -v C=$$(bawk '$$2=="Z85996.1" {print $$4}' $<) '$$3>=C' $^2 | gzip > $^3

greater_neg_in_common:
	matrix_reduce -t 'tpx_analysis/*/cCRE.tpx.best.complete.*-neg_pos.best_cut.gz' | bawk '$$6=="neg" {print $$3,$$2,$$4}' | cut -f 1 | symbol_count -d > $@
greater_neg: cCRE.bed
	matrix_reduce -t 'tpx_analysis/*/cCRE.tpx.best.complete.*-neg_pos.best_cut.gz' | bawk '$$6=="neg"' | translate -a <(echo -e "region\tpos_lncRNA"; cut -f 4,5 $<) 1 > $@ 
greater_neg_in_common.matrix: cCRE.bed
	matrix_reduce -t 'tpx_analysis/*/cCRE.tpx.best.complete.*-neg_pos.best_cut.gz' | bawk '$$6=="neg" {print $$3,$$2,$$4}' | tab2matrix -r region | translate -a <(echo -e "region\tpos_lncRNA"; cut -f 4,5 $<) 1 > $@

#tpx_analysis/%/cCRE.tpx.custom_t_pot: tpx_analysis/%/cCRE.tpx.gz cCRE.bed
#	bawk '$$Guanine_rate>=0.7 {print $$Duplex_ID,$$TTS_start,$$TTS_end}' $< | unhead | bedtools sort | bedtools merge | bawk '{print $$1,$$3-$$2}' | translate -a -r <(bawk '{print $$4,$$3-$$2}' $^2) 1 > $@

.META: cCRE.tpx.custom_t_pot.neg_pos_rand
	1	region	chirp_peak_1003;TERC
	2	covered_len	12
	3	total_len	804
	4	neg_pos	pos
	5	custom_t_pot	0.014925373134328358

%.tpx.summary.neg_pos: %.tpx.summary
	translate -a -k <(cut -f 4,6 $<) 1 < $^2 | bawk 'BEGIN{print "Duplex_ID","neg_pos","t_pot"} {print $$1,$$2,$$5}' > $@

%.tpx.custom_summary.pre: %.tpx cCRE.bed.fa TERC.fa
	unhead $< | cut -f 1,4 | symbol_count | translate -a -r <(fasta_length < $^2) 2 | translate -a -r <(fasta_length < $^3) 1 > $@

TERC-cCRE.bed.tpx.raw_%.custom_summary: TERC-cCRE.bed.tpx.raw_%.custom_summary.pre
	../../local/src/tiplexator_t_pot_norm.py -l 3 -L -1 < $< > $@
#TERC-cCRE.bed.tpx.summary.neg_pos: TERC.neg_pos_rand.bed TERC-cCRE.bed.tpx.raw.summary
#	translate -a -k <(cut -f 4,6 $<) 1 < $^2 | bawk 'BEGIN{print "Duplex_ID","neg_pos","t_pot"} {print $$1,$$2,$$5}' > $@

.META: TERC-cCRE.bed.tpx.raw_*.custom_summary
	1	ssRNA	TERC
	2	Duplex_ID	rand_chirp_peak_9;TERC
	3	txp_count	
	4	Duplex_lenght	780
	5	ssRNA_lenght	541
	6	nom_factor	174247920
	7	custom_t_pot	4.017264596329184e-08

TERC-cCRE.bed.tpx.raw_%.summary.custom_summary.neg_pos: TERC-cCRE.bed.tpx.raw_%.custom_summary.header_added TERC-cCRE.bed.tpx.raw_%summary.neg_pos
	translate -a -r <(cut -f 2,7 $<) 1 < $^2 > $@

TERC-cCRE.bed.tpx.raw%.covered_frac: TERC-cCRE.bed.tpx.raw% cCRE.bed
	bawk '{print $$4,$$5,$$6}' $< | unhead | bedtools sort | bedtools merge | bawk '{print $$1,$$3-$$2}' | stat_base -g -t | translate -a -r <(bawk '{print $$4,$$3-$$2}' $^2) 1 | bawk '{print $$1,$$2,$$2/$$3}' > $@

.META: TERC-cCRE.bed.tpx.raw*.covered_frac
	1	Duplex_ID
	2	TTS_covered_len
	3	TTS_covered_frac

TERC-cCRE.bed.tpx.raw%.neg_pos_rand: TERC.neg_pos_rand.bed TERC-cCRE.bed.tpx.raw%
	translate -a -r <(echo -e "Duplex_ID\tneg_pos"; cut -f 4,6 $<) 1 < $^2 > $@

TERC-cCRE.bed.tpx.raw_%.summary.clean: TERC-cCRE.bed.tpx.raw_%.summary
	bawk 'BEGIN{print "Duplex_ID","ssRNA","summary_total_tpx","summary_t_pot"} NR>1 {print $$1,$$2,$$3,$$4}' $< > $@

TERC-cCRE.bed.tpx.raw_%.summary.clean.covered_frac: TERC-cCRE.bed.tpx.raw_%.covered_frac.header_added TERC-cCRE.bed.tpx.raw_%.summary.clean
	translate -a $< 1 < $^2 > $@

TERC-cCRE.bed.tpx.raw_%.stability.best: TERC-cCRE.bed.tpx.raw_%.stability
	find_best 4 14 < $< > $@
TERC-cCRE.bed.tpx.raw_%.stability.tot_norm: TERC-cCRE.bed.tpx.raw_%.stability cCRE.bed
	cut -f 4,14 < $< | bsort | stat_base -g -t | translate -a -r <(bawk '{print $$4,$$3-$$2}' $^2) 1 | bawk '{print $$1,$$2,$$2/$$3}' > $@

TERC-cCRE.bed.tpx.raw_%.summary.clean.covered_frac.stability: TERC-cCRE.bed.tpx.raw_%.stability.best TERC-cCRE.bed.tpx.raw_%.summary.clean.covered_frac TERC-cCRE.bed.tpx.raw_%.stability.tot_norm
	translate -a -r <(bawk 'BEGIN{print "Duplex_ID","stability_best"} {print $$4,$$14}' $<) 1 < $^2 | translate -a -r <(bawk 'BEGIN{print "Duplex_ID","stability_tot","stability_norm"} {print}' $^3) 1 > $@

TERC-cCRE.bed.tpx.raw_%.summary.clean.covered_frac.stability.custom_t_pot.neg_pos_rand: TERC-cCRE.bed.tpx.raw_%.custom_summary.header_added TERC-cCRE.bed.tpx.raw_%.summary.clean.covered_frac.stability.neg_pos_rand
	translate -a -r -f 2 $< 1 < $^2 1 > $@





##################################
#
#	2D structure
#

%_lunp: %.fa
	conda activate rnaFoldRand_v0.1; \
	RNAplfold -W 200 -L 150 -u 10 -o < $<

#%_W200_L150_u10_0001_lunp: %.fa
#	conda activate rnaFoldRand_v0.1; \
#       RNAplfold -W 200 -L 150 -u 10 -o --id-prefix $*_W200_L150_u10 < $<

%_lunp.lastcol: %_lunp
	grep -v '#' $< | bawk '{print $$(NF)}' > $@

%_modif_zscore: %_lunp.lastcol
	conda activate pybigwig; \
        ../../local/src/distMedian.py < $< > $@


%_modif_zscore.bedGraph: %_modif_zscore
	grep -v '^#' $< | enumerate_rows -s 8 | bawk '{print "ssRNA",$$1-4,$$1+1-4,$$2}' > $@


%.single_nt_prob: %_lunp
	unhead -n 2 $< | cut -f1,2 > $@

%_modif_zscore.aboveMedian.bed: %_modif_zscore
	grep -v '^#' $< | bawk '$$1<0 {print "$*",NR-1+5,NR+5,$$1}' | bedtools merge > $@

%_ss20_unpaired_window.fa: RNAplfold/%_lunp.unpaired_window.modif_zscore %.fa
	PERC=$$(sort -n $< | awk '{all[NR] = $$0} END{print all[int(NR*0.2 - 0.5)]}'); fasta_mask <(bawk -v perc=$$PERC '$$1<perc {print "TERC",NR-1+4,NR+4}' $< | bedtools merge) < $^2 > $@
