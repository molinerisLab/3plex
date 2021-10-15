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

tpx_analysis.fisher_select_cutoff.ALLconditions:
	matrix_reduce -t 'tpx_analysis/*/cCRE.tpx.best.complete.*_neg_pos.fisher_select_cutoff' | tr ";" "\t" | cut -f1,2,4- > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.best: tpx_analysis.fisher_select_cutoff.ALLconditions
	find_best -r 1:2:3 10 < $< > $@

.META: tpx_analysis.fisher_select_cutoff.ALLconditions tpx_analysis.fisher_select_cutoff.ALLconditions.best
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

tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix: tpx_analysis.fisher_select_cutoff.ALLconditions.best
	bawk '$$3=="$*" {$$pvalue=$$pvalue+2.26261e-288; print $$lncRNA,$$condition";"$$cCRE,-log($$pvalue)/log(10)}' $< | tab2matrix -r GeneID > $@
tpx_analysis.fisher_select_cutoff.ALLconditions.matrix: tpx_analysis.fisher_select_cutoff.ALLconditions.best
	bawk '{$$pvalue=$$pvalue+2.26261e-288; print $$lncRNA,$$condition";"$$cCRE,-log($$pvalue)/log(10)}' $< | tab2matrix -r GeneID > $@

tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix.heatmap_euclid.pdf: tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix
	heatmap --pdf -e 13 -w 3 -n -c 12 -r 4 $< $@
tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix.heatmap_pearson.pdf: tpx_analysis.fisher_select_cutoff.ALLconditions.%.matrix
	heatmap --pdf -e 13 -w 3 -n -c 12 -r 4 -d correlation $< $@

GEP.count.exp_filter.ltmm.metadata.exp_genes_condition.matrix.selected_lncRNA.heatmap.pdf: GEP.count.exp_filter.ltmm.metadata.exp_genes_condition.matrix.selected_lncRNA.header_added
	heatmap -c 12 -e 13 -w 3 --pdf $< $^2
