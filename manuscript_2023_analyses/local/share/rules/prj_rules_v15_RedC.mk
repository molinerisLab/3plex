RADICLseq.ALL_ssRNA.bed.gz: ../../local/share/data/RADICLseq/GSE132190_mESC_NPM_ALL_significant.merge_distinct.matrix_reduce.merge_replicates.biotypes.gz
	zcat $< | bawk '$$5=="lincRNA" && $$6=="n1;n2" {print $$1,$$2,$$3,"RADICL_peak_"NR,$$4}' | gzip > $@

selected_ssRNA: RADICLseq.ALL_ssRNA.bed.gz
	zcat $< | cut -f5 | sort -u > $@

%.neg_pos.bed: %_onlypos.bed %_neg.bed
	cat $^ > $@

%_onlypos.bed: RADICLseq.ALL_ssRNA.bed.gz
	bawk '$$5=="$*" || $$5~/$*.intersec/ {print $$1,$$2,$$3,$$4";"$$5,$$5,"pos"}' $< > $@

%_neg.bed: RADICLseq.ALL_ssRNA.bed.gz %_onlypos.bed
	bawk '{print $$1,$$2,$$3,$$4";"$$5,"$*"}' $< > $@.tmp
	bedtools intersect -a $@.tmp -b $^2 -v \
	| bawk '{print $$0, "neg" }'> $@
	rm $@.tmp

%.neg_rand.excl.bed: $(GENCODE_DIR)/mm10.shuffle_blacklist.bed $(GENCODE_DIR)/gap.bed %_onlypos.bed
	cut -f -3 $^ |  bedtools sort | bedtools merge > $@

%.neg_rand.bed: %_onlypos.bed $(GENCODE_DIR)/chrom.info.no_alt %.neg_rand.excl.bed
	bedtools shuffle -i $< -g $^2 -excl $^3 | bawk '{$$4="rand_"$$4; $$6="neg"; print}' > $@

%.neg_pos_rand.bed: %_onlypos.bed %.neg_rand.bed
	cat $^ > $@




ALL_ssRNA=$(shell cat selected_ssRNA_id)

ALL_neg_pos_rand.bed: $(addsuffix .neg_pos_rand.bed,$(ALL_ssRNA))
	cut -f -4 $^ > $@
ALL_neg_pos.bed: $(addsuffix .neg_pos.bed,$(ALL_ssRNA))
	cut -f -4 $^ > $@