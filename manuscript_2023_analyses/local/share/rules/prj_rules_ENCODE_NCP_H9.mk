
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
hg38-ccREs.bed.fa.tpx.%_pos_neg.fisher_select_cutoff: hg38-ccREs.bed.fa.tpx.%_pos_neg
	bawk '{print $$1";"$$4,$$3,$$5}' $< | ./fisher_select_cutoff.py | tr ";" "\t" > $@

.META: hg38-ccREs.bed.fa.tpx.*_pos_neg.fisher_select_cutoff
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

hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fisher_select_cutoff: hg38-ccREs.bed.fa.tpx.NPC_H9_pos_neg.fisher_select_cutoff hg38-ccREs.bed.fa.tpx.H9_pos_neg.fisher_select_cutoff hg38-ccREs.bed.fa.tpx.NPCvsH9_pos_neg.fisher_select_cutoff
	matrix_reduce 'hg38-ccREs.bed.fa.tpx.*_pos_neg.fisher_select_cutoff' -l '$^' | fasta2tab | bawk '{print $$2~10,$$1}' > $@

.META: hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fisher_select_cutoff
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


hg38-cCRE_specific_H9-NPC: hg38-H9.neg_pos.bed hg38-NPC_H9.neg_pos.bed
	matrix_reduce 'hg38-*.neg_pos.bed' -l '$^' | fasta2tab | bawk '$$7=="pos" {print $$5,$$6,$$1}' | sed 's/NPC_H9/NPC/' | collapsesets 3 | grep -v ';' > $@
hg38-ccREs.bed.fa.tpx.best_score.specific_H9-NPC: hg38-cCRE_specific_H9-NPC hg38-ccREs.bed hg38-ccREs.bed.fa.tpx.best_score
	translate -a -r -k -d <(cat $< | translate -f 2 <(cut -f 4,5 $^2) 1) 2 < $^3 \
	| grep -v ';' > $@                      *perdiamo i ccRE che cambiano tipo ad H9 a NPC, ma sono pochi e nessun dELS*
hg38-ccREs.bed.fa.tpx.NPCvsH9_pos_neg.fisher_select_cutoff: hg38-ccREs.bed.fa.tpx.best_score.specific_H9-NPC
	bawk '{if($$5=="NPC"){$$5="pos"}else{$$5="neg"} print $$1";"$$4,$$3,$$5}' $< | fisher_select_cutoff  | tr ";" "\t" > $@

upregulated: /sto1/epigen/ENCODE_RNAseq/dataset/NPC_H9/DGE/selected_lncRNA
	bawk '{if($$6==1){up="NPC"}else{up="H9"} print $$7,up}' $< > $@
hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fisher_select_cutoff.upregulated: upregulated hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fisher_select_cutoff
	translate -a -r $< 1 < $^2 > $@

.META: hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fisher_select_cutoff.upregulated
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
	11	upregulated

hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fisher_select_cutoff.upregulated.best_pvalue: hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fisher_select_cutoff.upregulated ssRNA.len
	bawk 'NR>1 && $$cCRE_source=="NPCvsH9" {s=1; if($$oddsratio<1){s=-1}; print $$lncRNA";"$$cCRE";"$$upregulated, -log($$pvalue)/log(10), -s*log($$pvalue)/log(10)}' $< | find_best 1 2 | tr ";" "\t" \
	| translate -a -r $^2 1 > $@

.META: hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fisher_select_cutoff.upregulated.best_pvalue
	1	lncRNA
	2	cCRE
	3	upregulated
	4	score_abs
	5	score
	6	len

