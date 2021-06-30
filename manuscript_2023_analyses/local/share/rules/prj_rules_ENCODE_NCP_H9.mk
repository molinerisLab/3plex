
ssRNA.len: /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/gencode.v32.transcripts.fa.gz ../../local/share/data/selected_lncRNA.ENCODE_NPC_H9
	zcat $< | fasta2oneline | tr "|" "\t" | filter_1col 6 <(cut -f 1 $^2) | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,7 > $@

ssRNA.fa: $(BIOINFO_REFERENCE_ROOT)/gencode/dataset/$(GENCODE_SPECIES)/$(GENCODE_VERSION)/gencode.v32.transcripts.fa.gz ../../local/share/data/selected_lncRNA.ENCODE_NPC_H9
	zcat $< | fasta2oneline | tr "|" "\t" | filter_1col 6 <(cut -f 1 $^2) | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,10 | tab2fasta | fold > $@

#############################################################################
#
#	cCRE specific rules
#

hg38-ccREs.bed.fa.tpx.best_score: hg38-ccREs.bed.fa.tpx
	bawk 'NR>1 && $$Score>=10 {print $$sequence_id";"$$Duplex_ID,$$Score}' $< | find_best -s 1 2 | tr ";" "\t" > $@
hg38-ccREs.bed.fa.tpx.best_score.complete: hg38-ccREs.bed.fa.tpx.best_score hg38-ccREs.bed
	tab2matrix -e 0 -C <(cut -f 4 $^2) < $< | matrix2tab > $@
hg38-ccREs.bed.fa.tpx.%_pos_neg: hg38-ccREs.bed.fa.tpx.best_score.complete hg38-ccREs.bed hg38-%.neg_pos.bed
	translate <(cut -f 4,5 $^2) 2 < $< | translate -a -r -k <(cut -f 4- $^3) 2 > $@
hg38-ccREs.bed.fa.tpx.%_pos_neg.fischer_select_cutoff: hg38-ccREs.bed.fa.tpx.%_pos_neg
	bawk '{print $$1";"$$4,$$3,$$5}' $< | ./fischer_select_cutoff.py -a greater | tr ";" "\t" > $@

.META: hg38-ccREs.bed.fa.tpx.*_pos_neg.fischer_select_cutoff
	1	lncRNA
	2	cCRE
	3	score
	4	greater_positive
	5	greater_negative
	6	lower_positive
	7	lower_negative
	8	oddsratio
	9	pvalue

hg38-ccREs.bed.fa.tpx.ALL_pos_neg: hg38-ccREs.bed.fa.tpx.NPC_H9_pos_neg hg38-ccREs.bed.fa.tpx.H9_pos_neg
	matrix_reduce 'hg38-ccREs.bed.fa.tpx.*_pos_neg' -l '$^' | fasta2tab > $@

.META: hg38-ccREs.bed.fa.tpx.ALL_pos_neg
	1	cCRE_source	NPC_H9
	2	lncRNA		LINC00461
	3	cCRE_id		EH38E0065969
	4	score		10
	5	cCRE_type	dELS
	6	pos_neg		neg

hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fischer_select_cutoff: hg38-ccREs.bed.fa.tpx.NPC_H9_pos_neg.fischer_select_cutoff hg38-ccREs.bed.fa.tpx.H9_pos_neg.fischer_select_cutoff
	matrix_reduce 'hg38-ccREs.bed.fa.tpx.*_pos_neg.fischer_select_cutoff' -l '$^' | fasta2tab | bawk '{print $$2~10,$$1}' > $@

.META: hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fischer_select_cutoff
	1	lncRNA
	2	cCRE
	3	score
	4	greater_positive
	5	greater_negative
	6	lower_positive
	7	lower_negative
	8	oddsratio
	9	pvalue
	10	cCRE_source

