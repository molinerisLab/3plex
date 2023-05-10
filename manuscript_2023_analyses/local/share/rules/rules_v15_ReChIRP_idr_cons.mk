########################
#
#     General param
#

GENCODE_DIR=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/mmusculus/M25
NCPU?=28


###########################
# 
#      Docker param
#

DOCKER_GROUP?=$(shell stat i-c %g $(PWD))
DOCKER_DATA_DIR?=$(PRJ_ROOT)


FASTA_KTUP?=6

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




#################################
#
#	v8 ROC analysis
#

# matrix_reduce -t 'tpx_paramspace/*_*_*/*.neg_pos_rand.bed/*/*/*/*/*/*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.AUC_comp.gz'
tpx_paramspace_AUC_cmp.gz:
	matrix_reduce -t 'tpx_paramspace/*_*_*/*.neg_pos_rand.bed/*/*/*/*/*/*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.AUC_comp.gz' \
	| grep -v -w pred_1 | tr ";" "\t" \
	| perl -lane '$$,="\t"; @F=map{s/.*\~//; $$_} @F; print @F' \
	| cut -f 1-3,5-  \
	| bawk '{print $$0; print $$1~9,$$11,$$10,$$13,$$12,$$14}' | sort | uniq | gzip > $@

tpx_paramspace_AUC.gz: tpx_paramspace_AUC_cmp.gz	
	zcat $< | cut -f -10,12 | uniq | gzip > $@

.META: tpx_paramspace_AUC_cmp.gz
	1	ssRNA	AC018781.1
	2	single_stranddnes_cutoff	ss0
	3	RNAplfold_window	singleNt
	4	min_length
	5	max_length
	6	error_rate
	7	guanine_rate
	8	filter_repeat
	9	consecutive_errors
	10	predictor1
	11	predictor2
	12	AUC1
	13	AUC2
	14	pvalue

.META: tpx_paramspace_AUC.gz
	1	ssRNA	AC018781.1
	2	single_stranddnes_cutoff	ss0
	3	RNAplfold_window	singleNt
	4	min_length
	5	max_length
	6	error_rate
	7	guanine_rate
	8	filter_repeat
	9	consecutive_errors
	10	predictor
	11	AUC
	12	Pvalue	Mann Withey
	13	PvalueAdj BH

tpx_paramspace_AUC.all_human.gz: ../v8_ChIRP_neg_rand/tpx_paramspace_AUC.gz ../v8.2_ReChIRP_idr_cons/tpx_paramspace_AUC.gz ../v8.3_ReChIRP_overlap/tpx_paramspace_AUC.gz ../v8.6_ReChIRP_idr_overlap_top1000/tpx_paramspace_AUC.gz
	matrix_reduce -t -l '$^' '../*/tpx_paramspace_AUC.gz' | gzip > $@

.META: tpx_paramspace_AUC.all_human.gz
	1	selected_peaks
	2	ssRNA
	3	single_stranddnes_cutoff
	4	RNAplfold_window
	5	min_length
	6	max_length
	7	error_rate
	8	guanine_rate
	9	filter_repeat
	10	consecutive_errors
	11	predictor
	12	AUC

PROB__fitted_model_evaluation_fixed_param: tpx_paramspace_AUC_cmp.gz
	echo -e "lncRNA\tCSS AUC\tT_POT AUC\tP-value" > $@
	bawk '$$RNAplfold_window=="singleNt" && $$guanine_rate==40 && $$error_rate==20 && $$max_length==-1 && $$min_length==10 && $$single_stranddnes_cutoff=="ss0" && $$filter_repeat=="off" && $$consecutive_errors==3 && $$predictor1=="PROB__fitted_model" && $$predictor2=="t_pot_norm" {pn=1; best_AUC=$$AUC1; if($$AUC1<$$AUC2){pn=2; best_AUC=$$AUC2} print $$0,pn,best_AUC}' $< \
	| round_table -p 3 | sed 's/PROB__fitted_model/CSS/; s/t_pot_norm/T_POT/' | cut -f 1,12,13,14 >> $@

PROB__fitted_model_evaluation_best_param: tpx_paramspace_AUC_cmp.gz
	echo -e "lncRNA\tCSS AUC\tT_POT AUC\tP-value" > $@
	bawk '$$predictor1=="PROB__fitted_model" && $$predictor2=="t_pot_norm" {pn=1; best_AUC=$$AUC1; if($$AUC1<$$AUC2){pn=2; best_AUC=$$AUC2} print $$0,pn,best_AUC}' $< \
	| find_best 1 16 \
	| round_table -p 3 | sed 's/PROB__fitted_model/CSS/; s/t_pot_norm/T_POT/' | cut -f 1,12,13,14,15 >> $@

PROB__fitted_model_unpairedWindow_evaluation_best_param: tpx_paramspace_AUC_cmp.gz
	bawk '$$RNAplfold_window=="unpairedWindow" && $$predictor1=="PROB__fitted_model"' $< | find_best 1 12 > $@

bestAUC_params.tsv: tpx_paramspace_AUC_cmp.gz
	bawk '$$RNAplfold_window=="unpairedWindow" && $$predictor1=="PROB__fitted_model"' $< | find_best 1 12 > $@

#	bawk 'BEGIN{print "GeneID","RNA_ss_cutoff","RNAplfold_window","min_length","max_length","error_rate","guanine_rate","filter_repeat","consecutive_errors"} {print}' > $@

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.bestAUC_params.matrix_reduce: bestAUC_params.tsv
	matrix_reduce -t 'tpx_paramspace/*_*_unpairedWindow/*.neg_pos_rand.bed/min_length~10/max_length~*/error_rate~20/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic' | filter_1col 1 <(bawk '{print $$1";"$$2";"$$1";"$$5";"$$7";"$$8";"$$9}' $< | unhead) | tr ";" "\t" | bawk '{print $$1,$$8";"$$9,$$10~25}' | grep -v "Stability" > $@

.META: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.bestAUC_params.matrix_reduce
	1	ssRNA	AC018781.1
	2	DuplexID	chirp_peak_31;AC018781.1
	3       tpx_count_custom
	4	neg_pos
	5	tpx_count_standard
	6	t_pot_norm	
	7	Duplex_length
	8	oligolength
	9	triplexator_norm_factor
	10	t_pot_custom
	11	TTS_covered_len
	12	TTS_covered_frac
	13	Stability_best
	14	Stability_tot_overcount
	15	Stability_tot_undercount
	16	Stability_norm_overcount
	17	Stability_norm_undercount
	18	PROB__fitted_model


parameter_evaluation-max_length: tpx_paramspace_AUC_cmp.gz
	zcat $< | round_table -p 3 | find_best -m 1 12 | cut -f -10,12 | sort | uniq | collapsesets 5 | collapsesets 10 > $@

bestParams_bestPredictor.tsv: tpx_paramspace_AUC_cmp.gz
	bawk '$$5==-1' $< | find_best -m 1 12 | cut -f -10,12 | sort | uniq | collapsesets 3 | bawk '$$1!="LINC01605" {split($$3,a,";"); print $$1";"$$2";"a[1]";"$$1";"$$4";"$$7";"$$8";"$$9,$$1,$$10}' > $@

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.bestParams_bestPredictor.matrix_reduce: bestParams_bestPredictor.tsv
	matrix_reduce -t 'tpx_paramspace/*_*_*/*.neg_pos_rand.bed/min_length~*/max_length~-1/error_rate~20/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic' | filter_1col 1 <(cut -f1 $<) | tr ";" "\t" | bawk '{print $$1,$$9";"$$10,$$11~26}' | grep -v "Stability" > $@



selected_ssRNA.conditions.clean: selected_ssRNA.conditions
	filter_1col 1 <(cut -f 1 $< | symbol_count | bawk '$$2==1 {print $$1}') < $< > $@
selected_ssRNA.conditions.up_down.balanced.clean: selected_ssRNA.conditions.up_down.balanced
	filter_1col -v 1 <(bawk '$$4==1 {print $$1}' $< | symbol_count | bawk '$$2>1 {print $$1}') < $< > $@

cCRE-%.type: cCRE-%.matrix
	tr ";" "\t" < $< | unhead | bawk '{type=$$5; if(type=="Low-DNase"){type=$$6} print $$4,type}' > $@

cCRE-%.type.bed: cCRE-%-vs-H9.type cCRE-%.bed
	translate -a $< 4 < $^2 > $@
cCRE-%.neg_pos.bed: cCRE-%-vs-H9.matrix cCRE-%.type.bed
	translate -a -r <(tr ";" "\t" < $< | unhead | cut -f 4,5 | bawk '$$2!="Low-DNase" {print $$1,"pos"} $$2=="Low-DNase" {print $$1,"neg"}') 4 < $^2 > $@

cCRE-%.type.bed: cCRE-%-vs-H1.type cCRE-%.bed
	translate -a $< 4 < $^2 > $@
cCRE-%.neg_pos.bed: cCRE-%-vs-H1.matrix cCRE-%.type.bed
	translate -a -r <(tr ";" "\t" < $< | unhead | cut -f 4,5 | bawk '$$2!="Low-DNase" {print $$1,"pos"} $$2=="Low-DNase" {print $$1,"neg"}') 4 < $^2 > $@

tpx_paramspace.fisher_select_cutoff.matrix:
	matrix_reduce -t 'tpx_paramspace/*_ss*_unpairedWindow/cCRE-*.bed/min_length~*/max_length~*/error_rate~*/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.*.type.neg_pos.fisher_select_cutoff' | tr ";" "\t" > $@



tpx_paramspace.fisher_select_cutoff.matrix.best: tpx_paramspace.fisher_select_cutoff.matrix
	find_best -r 1:2:3:4:5:6:7:8:9:10:11 18 < $< | tr ";" "\t" > $@

.META: tpx_paramspace.fisher_select_cutoff.matrix tpx_paramspace.fisher_select_cutoff.matrix.best
	1	ssRNA	AC004797.1
	2	secondary_structure	0
	3	condition	HEP_H9
	4	min_length	10
	5	max_length	-1
	6	error_rate	20
	7	guanine_rate	10
	8	filter_repeat	off
	9	consecutive_errors	1
	10	score_type	stability.norm_undercount
	11	cCRE_type	dELS
	12	score_cutoff	0.98
	13	greater_positive
	14	lower_positive
	15	greater_negative
	16	lower_negative
	17	oddsratio	2.345776
	18	Pvalue	2.07884e-24
	19	TPXcCRE_score	23.682180

tpx_paramspace.fisher_select_cutoff.matrix.best.up_down: tpx_paramspace.fisher_select_cutoff.matrix.best.header_added selected_ssRNA.conditions.up_down.balanced.clean
	bawk '{print $$ssRNA";"$$condition,$$2,$$4~19}' $< | translate -a -k <(bawk 'BEGIN{print "ssRNA;condition", "staminal", "mark_seqc", "logFC"} {print $$1";"$$2,$$3~5}' $^2) 1 | tr ";" "\t" > $@

%.fa.nin: %.fa
	makeblastdb -in $< -dbtype nucl

%.fa.hoogsteen: %.fa
	../../local/src/word_count_RNA_hoogsteen.pl < $< > $@


AC009055.2.fa.hoogsteen-cCRE-ALLcond.fasta_out: AC009055.2.fa.hoogsteen cCRE-ALLcond.bed.fa
	conda activate /sto1/ref/miniconda2/envs/bit_rnaseq_2.8/envs/blast/;\	
	fasta36 -n -b 1000000 -d 10000000 -z -1 -m 9C -T $(NCPU) $< $^2 $(FASTA_KTUP) > $@

	#fasta36 -3 -n -b 1000000 -d 10000000 -z -1 -m 2 -T $(NCPU) $< $^2 > $@

%.fa.hoogsteen-EH38E1935892.fasta_m10: %.fa.hoogsteen EH38E1935892.fa
	conda activate /sto1/ref/miniconda2/envs/bit_rnaseq_2.8/envs/blast/;\
	fasta36 -n -b 1000000 -d 10000000 -T 28 -m 10 $< $^2 3 > $@
%.fa.hoogsteen-EH38E1935892.fasta_m9c: %.fa.hoogsteen EH38E1935892.fa
	conda activate /sto1/ref/miniconda2/envs/bit_rnaseq_2.8/envs/blast/;\
	fasta36 -n -b 1000000 -d 10000000 -T 28 -m 9c $< $^2 3 > $@

%.fa.hoogsteen-EH38E1935892.fasta_m10.parsed: %.fa.hoogsteen-EH38E1935892.fasta_m10
	conda activate /sto1/ref/miniconda2/envs/bit_rnaseq_2.8/envs/blast/;\
	../../local/src/fasta2tpx.py < $< > $@
%.fa.hoogsteen-EH38E1935892.fasta_m10.parsed.tpx_pre: %.fa.hoogsteen-EH38E1935892.fasta_m10.parsed
	../../local/src/break_aln.pl -l 1 < $< > $@

AC003681.1_EH38E1935892_triplexator.tpx: EH38E1935892.fa AC003681.1.fa
	docker run -u `id -u`:10001 --rm -v /sto1:/sto1 -v /sto1:/sto1 triplexator:v1.3.2_l6 bash -c "cd /sto1/epigen/TPXcCRE/dataset/v14_all_cCRE_conditions; \
	triplexator -l 8 -L -1 -e 20 -g 10 -fr off -c 1 -fm 0 -of 1 -o $@ -rm 2 -p 4 -ss $^2 -ds $<"

CALML3-AS1_EH38E1310212_triplexator.tpx: EH38E1310212.fa CALML3-AS1.fa
	docker run -u `id -u`:10001 --rm -v /sto1:/sto1 -v /sto1:/sto1 triplexator:v1.3.2_l6 bash -c "cd /sto1/epigen/TPXcCRE/dataset/v14_all_cCRE_conditions; \
	triplexator -l 8 -L -1 -e 20 -g 10 -fr off -c 1 -fm 0 -of 1 -o $@ -rm 2 -p 4 -ss $^2 -ds $<"

tpx_paramspace_AUC.all_human.matrix: tpx_paramspace_AUC.all_human.gz
	bawk '{print $$1";"$$2,$$3~12}' $< | grep -v MIR503HG | find_best 1 11 | tr ";" "\t" | cut -f -2,12 | tab2matrix -t > $@




#############################
#
# tpx_paramspace analysis
#

tpx_paramspace_AUC.all_%_noSingleNt.gz: ../raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.AUC_comp.ALL_versions_peralta.gz %_lncRNA
	cp $< $@.tmp.gz; \
	zcat $@.tmp.gz | filter_1col 2 $^2 | cut -f-9,11 | sort | uniq | gzip > $@; \
	zgrep t_pot_norm $@.tmp.gz | filter_1col 2 $^2 | cut -f-8,10,12 | sort | uniq | gzip >> $@; \
	rm $@.tmp.gz

tpx_paramspace_AUC.all_%_noSingleNt.method.Npeaks.gz: tpx_paramspace_AUC.all_%_noSingleNt.gz ../../local/share/data/all_reproducibility.all_bed.regionPeak.counts_matrix
	zcat $< | grep -v LINC01605 | grep -v HOTTIP | translate -a <(cut -f2-5,7,10 $^2) 2 | gzip > $@

tpx_paramspace_AUC.all_human_noSingleNt.method.Npeaks_version.gz: tpx_paramspace_AUC.all_human_noSingleNt.method.Npeaks.gz
	bawk '{if($$1 == "v8.2_ReChIRP_idr_cons"){print $$1~2,$$4,$$8~15}else if($$1 == "v8.3_ReChIRP_overlap"){print $$1~2,$$6,$$8~15} else if($$1 == "v8.6_ReChIRP_idr_overlap_top1000"){print $$1~2,1000,$$8~15} else if($$1 == "v8_ChIRP_neg_rand") {print $$1~2,$$7~15} else if($$1 == "v8.8_ReChIRP_idr_opt") {print $$1~2,$$5,$$8~15}}' $< | gzip > $@

tpx_paramspace_AUC.all_mouse_noSingleNt.method.Npeaks_version.gz: tpx_paramspace_AUC.all_mouse_noSingleNt.method.Npeaks.gz
	bawk '{if($$1 == "v8.4_ReChIRP_idr_mouse"){print $$1~2,$$4,$$8~15}else if($$1 == "v8.5_ReChIRP_overlap_mouse"){print $$1~2,$$6,$$8~15} else if($$1 == "v8.7_ReChIRP_idr_overlap_top1000_mouse"){print $$1~2,1000,$$8~15} else if($$1 == "v8.9_ReChIRP_idr_opt_mouse") {print $$1~2,$$5,$$8~15} else if($$1 == "v8.11_ReChIRP_bed_mouse"){print $$1~2,$$7~15}}' $< | gzip > $@

.META: tpx_paramspace_AUC.all_*_noSingleNt.method.Npeaks_version.gz
	1	peak_method	v8.2_ReChIRP_idr_cons
	2	ssRNA	7SK
	3	n_peaks	249
	4	singleStrandedness	ss0
	5	min_len	10
	6	error_rate	10
	7	guanine_rate	10
	8	repeat_filter	off
	9	consecutive_error	1
	10	predictor	PROB__fitted_model
	11	AUC	0.5

%.anova_exclude.gz: ./anova_perm.sh
	export OPENBLAS_NUM_THREADS=1; seq 1 10000 | parallel -j 50 --group $< > $@.tmp; \
	grep -v '^>' $@.tmp | gzip > $@; \
	rm $@.tmp

FORMULA=AUC~peak_method+ssRNA+n_peaks+singleStrandedness+min_len+error_rate+guanine_rate+repeat_filter+consecutive_error+predictor
%.true_prediction: %.gz
	bmodel -i -a $< $(FORMULA) | grep -v '^>' > $@









##############################
#
# Da revisionare
#

tpx_paramspace_AUC.all_%_noSingleNt.method.Npeaks_version.qc_params.gz: tpx_paramspace_AUC.all_%_noSingleNt.method.Npeaks_version.gz /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/qc_params.tsv
	zcat $< | translate -a -k <(cut -f2- $^2) 2 | gzip > $@
#	zcat $< | translate -a -k <(cut -f2- $^2 | grep -v 'AC018781') 2 | gzip > $@
ssRNA_frip_type.%.tsv: tpx_paramspace_AUC.all_%_noSingleNt.method.Npeaks_version.qc_params.gz
	zcat $< | cut -f2,4,5,8,9 | sort | uniq | bawk '{a = $$1; print a";overlap."$$2,a";overlap."$$3,a";idr."$$4,a";idr."$$5}' | tr "\t" "\n" > $@
tpx_paramspace_AUC.all_%_noSingleNt.method.Npeaks_version.qc_params.frip.gz: tpx_paramspace_AUC.all_%_noSingleNt.method.Npeaks_version.qc_params.gz /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/ssRNA_firp_samples.tsv ssRNA_frip_type.%.tsv
	zcat $< | translate -a -d <(translate -a -v -e NA $^2 1 < $^3 | tr ";" "\t") 2 | tr ";" "\t" | cut -f1-11,14,15,18- | gzip > $@

.META: tpx_paramspace_AUC.all_*_noSingleNt.method.Npeaks_version.qc_params.frip.gz
	1	peak_method	v8.2_ReChIRP_idr_cons
	2	ssRNA	7SK
	3	overlap.opt_set	overlap.pooled-pr1_vs_pooled-pr2
	4	overlap.opt.frip	0.3111992763582027
	5	overlap.cons_set	overlap.rep1_vs_rep2
	6	overlap.cons.frip	0.033778479652571976
	7	idr.opt_set	idr.pooled-pr1_vs_pooled-pr2
	8	idr.opt.frip	0.22089077906436314
	9	idr.cons_set	idr.rep1_vs_rep2
	10	idr.cons.frip	0.004493715610772667
	11	paired_end	False
	12	overlap.rescue_ratio	11.979719917012448
	13	overlap.self_consistency_ratio	1.3578253351920213
	14	idr.rescue_ratio	482.2570281124498
	15	idr.self_consistency_ratio	1.376415139595126
	16	exp_method	ChIRP-seq
	17	n_peaks	249
	18	singleStrandedness	ss0
	19	min_len	10
	20	max_len	-1
	21	error_rate	10
	22	guanine_rate	10
	23	repeat_filter	off
	24	consecutive_error	1
	25	predictor	PROB__fitted_model
	26	AUC	0.5

tpx_paramspace_AUC.all_human_noSingleNt.method.Npeaks_version.qc_params.frip_version.gz: tpx_paramspace_AUC.all_human_noSingleNt.method.Npeaks_version.qc_params.frip.gz
	bawk '{if($$1 == "v8.2_ReChIRP_idr_cons"){print $$1~2,$$10,$$14~26}else if($$1 == "v8.3_ReChIRP_overlap"){print $$1~2,$$6,$$12,$$13,$$16~26} else if($$1 == "v8.6_ReChIRP_idr_overlap_top1000"){print $$1~2,($$4+$$6+$$8+$$10)/4,($$12+$$14)/2,($$13+$$15)/2,$$16~26} else if($$1 == "v8_ChIRP_neg_rand") {print $$1~2,"NA","NA","NA",$$16~26} else if($$1 == "v8.8_ReChIRP_idr_opt"){print $$1~2,$$8,$$14~26}}' $< | expandsets 3 4 5 | gzip > $@

.META: tpx_paramspace_AUC.all_*_noSingleNt.method.Npeaks_version.qc_params.frip_version.gz
	1	peak_method	v8.6_ReChIRP_idr_overlap_top1000
	2	ssRNA	7SK
	3	frip	0.3111992763582027
	4	rescue_ratio	11.979719917012448
	5	self_consistency_ratio	1.3578253351920213
	6	exp_method	ChIRP-seq
	7	n_peaks	1000
	8	singleStrandedness	ss0
	9	min_len	10
	10	max_len	-1
	11	error_rate	10
	12	guanine_rate	10
	13	repeat_filter	off
	14	consecutive_error	1
	15	predictor	PROB__fitted_model
	16	AUC	0.5

tpx_paramspace_AUC.all_%_noSingleNt.method.Npeaks_version.qc_params.frip_version.finale.gz: tpx_paramspace_AUC.all_%_noSingleNt.method.Npeaks_version.qc_params.frip_version.gz
	zcat $< | cut -f1-5,7-9,11- | gzip > $@

.META: tpx_paramspace_AUC.all_*_noSingleNt.method.Npeaks_version.qc_params.frip_version.finale.gz
	1	peak_method	v8.2_ReChIRP_idr_cons
	2	ssRNA	7SK
	3	frip	0.004493715610772667
	4	rescue_ratio	482.2570281124498
	5	self_consistency_ratio	1.376415139595126
	6	n_peaks	249
	7	singleStrandedness	ss0
	8	min_len	10
	9	error_rate	10
	10	guanine_rate	10
	11	repeat_filter	off
	12	consecutive_error	1
	13	predictor	PROB__fitted_model
	14	AUC	0.5

tpx_paramspace_AUC.all_human_noSingleNt.method.Npeaks_version.qc_params.frip_version.finale.true_prediction: tpx_paramspace_AUC.all_human_noSingleNt.method.Npeaks_version.qc_params.frip_version.finale.header_added.gz
	bmodel -i -a $< 'AUC~peak_method+ssRNA+frip+rescue_ratio+self_consistency_ratio+n_peaks+singleStrandedness+min_len+error_rate+guanine_rate+repeat_filter+consecutive_error+predictor' | grep -v '^>' > $@

tpx_paramspace_AUC.idr_overlap_top1000.best_general_params.gz:
	matrix_reduce -t 'tpx_paramspace/*_ss0_unpairedWindow/*.neg_pos_rand.bed/min_length~8/max_length~-1/error_rate~20/guanine_rate~40/filter_repeat~off/consecutive_errors~1/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.gz' | grep -v Duplex_ID | bawk '{split($$1,a,";"); print a[1],$$2~18}' | gzip > $@

#tpx_paramspace_AUC.all_human_noSingleNt.gz: ../v8.2_ReChIRP_idr_cons/tpx_paramspace_AUC.gz ../v8.3_ReChIRP_overlap/tpx_paramspace_AUC.gz ../v8.6_ReChIRP_idr_overlap_top1000/tpx_paramspace_AUC.gz ../v8.8_ReChIRP_idr_opt/tpx_paramspace_AUC.gz
#	matrix_reduce -t -l '$^' '../*/tpx_paramspace_AUC.gz' | bawk '$$4 != "singleNt" {print $$1~3,$$5~13}' | gzip > $@
#tpx_paramspace_AUC.all_mouse_noSingleNt.gz: ../v8.4_ReChIRP_idr_mouse/tpx_paramspace_AUC.gz ../v8.5_ReChIRP_overlap_mouse/tpx_paramspace_AUC.gz ../v8.7_ReChIRP_idr_overlap_top1000_mouse/tpx_paramspace_AUC.gz ../v8.9_ReChIRP_idr_opt_mouse/tpx_paramspace_AUC.gz
#	matrix_reduce -t -l '$^' '../*/tpx_paramspace_AUC.gz' | bawk '$$4 != "singleNt" {print $$1~3,$$5~13}' | gzip > $@
#tpx_paramspace_AUC.idr_overlap_top1000.best_general_params.gz:
#	matrix_reduce -t 'tpx_paramspace/*_ss*_unpairedWindow/*.neg_pos_rand.bed/min_length~*/max_length~-1/error_rate~*/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.gz' | grep -v Duplex_ID | bawk '{split($$1,a,";"); print a[1],$$2~18}' | gzip > $@

.META: tpx_paramspace_AUC.idr_overlap_top1000.best_general_params.gz
	1	ssRNA
	2	Duplex_ID
	3	tpx_count_custom
	4	neg_pos
	5	tpx_count_standard
	6	t_pot_norm
	7	Duplex_length
	8	oligolength
	9	triplexator_norm_factor
	10	t_pot_custom
	11	TTS_covered_len
	12	TTS_covered_frac
	13	Stability_best
	14	Stability_tot_overcount
	15	Stability_tot_undercount
	16	Stability_norm_overcount
	17	Stability_norm_undercount
	18	PROB__fitted_model


#best_single_params.matrix.gz: best_parameter_settings.selected_files.tsv
#	@echo 'The -l argument does not work, I pasted the lines of the file manually'
#	matrix_reduce -t -l '<(cat $< | tr "\n" " ")' '../*/tpx_paramspace/*_*_unpairedWindow/*.neg_pos_rand.bed/min_length~*/max_^Cngth~-1/error_rate~*/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.gz' | bawk '{split($$1,a,";"); print a[2], $$2~18}' | grep -v 'Duplex_ID' | gzip > $@

#best_single_params.peak_method.matrix.gz:
#	@echo 'The -l argument does not work, I pasted the lines of the file manually'
#	matrix_reduce -t -l '<(ls $$(bawk '{print "../*/tpx_paramspace/"$$2"_"$$3"_unpairedWindow/"$$2".neg_pos_rand.bed/min_length~"$$4"/max_length~-1/error_rate~"$$6"/guanine_rate~"$$7"/filter_repeat~"$$8"/consecutive_errors~"$$9"/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.AUC_comp.gz"}' best_parameter_settings.tsv) | tr "\n" " ")' '../*/tpx_paramspace/*_*_unpairedWindow/*.neg_pos_rand.bed/min_length~*/max_length~-1/error_rate~*/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.AUC_comp.gz' | bawk '{split($$1,a,";"); print a[2], $$2~18}' | grep -v AUC | gzip > $@
#best_single_params.predictors.matrix.gz:
#	@echo 'The -l argument does not work, I pasted the lines of the file manually'
#	matrix_reduce -t -l '<(ls $$(bawk '{print "../*/tpx_paramspace/"$$2"_"$$3"_unpairedWindow/"$$2".neg_pos_rand.bed/min_length~"$$4"/max_length~-1/error_rate~"$$6"/guanine_rate~"$$7"/filter_repeat~"$$8"/consecutive_errors~"$$9"/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.AUC_comp.gz"}' best_parameter_settings.tsv) | tr "\n" " ")' '../*/tpx_paramspace/*_*_unpairedWindow/*.neg_pos_rand.bed/min_length~*/max_length~-1/error_rate~*/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.AUC_comp.gz' | bawk '{split($$1,a,";"); print a[2], $$2~18}' | grep -v AUC | gzip > $@


#best_parameter_settings.tsv: tpx_paramspace_AUC.all_human_noSingleNt.gz
#	zcat $< | grep -v LINC01605 | grep -v AC018781 | find_best 2 11 > $@
best_parameter_settings.selected_files.tsv: best_parameter_settings.tsv
	bawk '{print "../"$$1"/tpx_paramspace/"$$2"_"$$3"_unpairedWindow/"$$2".neg_pos_rand.bed/min_length~"$$4"/max_length~-1/error_rate~"$$6"/guanine_rate~"$$7"/filter_repeat~"$$8"/consecutive_errors~"$$9"/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.gz"}' $< > $@
best_single_params.matrix.gz:
	@echo 'Manually added v8_ChIRP_neg_rand values'


tpx_paramspace_AUC.all_mouse_noSingleNt.method.Npeaks_version.qc_params.frip_version.gz: tpx_paramspace_AUC.all_mouse_noSingleNt.method.Npeaks_version.qc_params.frip.gz
	bawk '{if($$1 == "v8.4_ReChIRP_idr_mouse"){print $$1~2,$$10,$$14~26}else if($$1 == "v8.5_ReChIRP_overlap_mouse"){print $$1~2,$$6,$$12,$$13,$$16~26} else if($$1 == "v8.7_ReChIRP_idr_overlap_top1000_mouse"){print $$1~2,($$4+$$6+$$8+$$10)/4,($$12+$$14)/2,($$13+$$15)/2,$$16~26} else if($$1 == "v8.8_ReChIRP_idr_opt"){print $$1~2,$$8,$$14~26}}' $< | expandsets 3 4 5 | gzip > $@

tpx_paramspace_AUC.all_mouse_noSingleNt.method.Npeaks_version.gz: tpx_paramspace_AUC.all_mouse_noSingleNt.method.Npeaks.gz
	bawk '{if($$1 == "v8.4_ReChIRP_idr_mouse"){print $$1~2,$$4,$$8~15}else if($$1 == "v8.5_ReChIRP_overlap_mouse"){print $$1~2,$$6,$$8~15} else if($$1 == "v8.7_ReChIRP_idr_overlap_top1000_mouse"){print $$1~2,1000,$$8~15} else if($$1 == "v8.9_ReChIRP_idr_opt_mouse") {print $$1~2,$$5,$$8~15} else if($$1 == "v8.11_ReChIRP_bed_mouse"){print $$1~2,$$7~15}}' $< | gzip > $@

tpx_paramspace_AUC.all_human_noSingleNt_noBed.gz:
	@echo 'Equivalent of tpx_paramspace_AUC.all_human_noSingleNt.gz on peralta'
best_parameter_settings.tsv: tpx_paramspace_AUC.all_human_noSingleNt.gz
	zcat $< | grep -v LINC01605 | grep -v AC018781 | find_best 2 11 > $@

# best_parameter_settings.selected_files.tsv: best_parameter_settings.tsv
# 	bawk '{print "../"$$1"/tpx_paramspace/"$$2"_"$$3"_unpairedWindow/"$$2".neg_pos_rand.bed/min_length~"$$4"/max_length~-1/error_rate~"$$6"/guanine_rate~"$$7"/filter_repeat~"$$8"/consecutive_errors~"$$9"/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.gz"}' $< > $@

best_single_params.matrix.gz:
	@echo 'Manually added v8_ChIRP_neg_rand values'

#tpx_paramspace_AUC.%.gz:
#	matrix_reduce -t 'tpx_paramspace/*_ss*_unpairedWindow/*.neg_pos_rand.bed/min_length~*/max_length~-1/error_rate~*/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.gz' | grep -v Duplex_ID | tr ";" "\t" | bawk '{print $$1,$$2,$$4~9";"$$10~26}' | gzip > $@

#tpx_paramspace_AUC.gz:
#	matrix_reduce -t 'tpx_paramspace/*_ss*_unpairedWindow/*.neg_pos_rand.bed/min_length~*/max_length~*/error_rate~*/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.AUC' | tr ";" "\t" | pvalue_correct -a -c 12 | gzip > $@

ssRNA.neg_pos_rand.bed:
	zgrep 'v8.2_ReChIRP_idr_cons' ../../local/share/data/ALL_v8.neg_pos_rand.bed.gz | bawk '{print $3~6";"$7,$8,$9 > $2".neg_pos_rand.bed"}'
raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.0_8_1_20_70_off_3.gz:
	matrix_reduce -t 'tpx_paramspace/*_ss0_unpairedWindow/*.neg_pos_rand.bed/min_length~8/max_length~-1/error_rate~20/guanine_rate~70/filter_repeat~off/consecutive_errors~3/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.gz' | \
	bawk 'NR==1 || $$2!="Duplex_ID" {split($$1,a,";"); print a[1],$$0}' | cut -f1,3- | gzip > $@
raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.0_8_1_20_70_off_3.AUC_cmp.tsv: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.0_8_1_20_70_off_3.gz
	zcat $< | ../../local/src/ROC.R neg_pos Stability_norm_undercount t_pot_norm Stability_best Score_best -d "<" -O raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.0_8_1_20_70_off_3.roc > $@
raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.10_10_1_20_10_off_3.gz:
	matrix_reduce -t 'tpx_paramspace/*_ss10_unpairedWindow/*.neg_pos_rand.bed/min_length~10/max_length~-1/error_rate~20/guanine_rate~10/filter_repeat~off/consecutive_errors~3/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.gz' | bawk 'NR==1 || $$2!="Duplex_ID" {split($$1,a,";"); print a[1],$$0}' | cut -f1,3- | gzip > $@
method_cmp.matrix.AUC_cmp.tsv: method_cmp.matrix.gz
	$(CONDA_ACTIVATE) pROC_Env;\
	zcat $< | ../../local/src/ROC.R neg_pos $$(zcat $< | head -n1 | cut -f3- | tr "\t" " ") -d "<" > $@

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.matrix.gz: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.gz selected_ssRNA.triplex_ssRNA
	bawk 'NR>1{split($$1,a,";"); id="d3plex_"a[4]"_"a[6]"_"a[7]"_"a[8]"_"a[9]"_"a[2]; print a[1],$$4";"$$2,id"_Stability_best",$$13; print a[1],$$4";"$$2,id"_Stability_norm",$$18}' $< | filter_1col 1 $^2 | cut -f2- | tab2matrix -r "pos_neg;peak;ssRNA" | bawk '{split($$1,a,";"); print a[1],$$0}' | cut -f 1,3- | gzip > $@
raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.matrix.auc_no_cmp.gz: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.matrix.gz ../../local/src/ROC_no_comp.R
	set +u; source /opt/conda/miniconda3/etc/profile.d/conda.sh; conda activate ; conda activate  pROC_Env;\
	zcat $< | $^2 pos_neg $$(zcat $< | head -n1 | cut -f2- | tr "\t" " ") -d "<" | gzip > $@
raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.matrix.auc_no_cmp.method.gz: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.matrix.auc_no_cmp.gz
	bawk 'BEGIN{print"shuffling","tpx_method","auc","species"}{print "genomic_regions",$$1,$$2,"human"}' $< | gzip > $@
