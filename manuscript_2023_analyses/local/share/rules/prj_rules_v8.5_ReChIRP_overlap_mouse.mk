ChIRP.bed.split.gz: ../../local/share/data/ReChIRP/mmusculus/overlap/overlap.conservative_peak.regionPeak.qTh005.mouse_matrix.gz
	zcat $< | bedtools sort | gzip > $@

%.neg_pos.bed: %_onlypos.bed %_neg.bed
	cat $^ > $@

%_onlypos.bed: ChIRP.bed.split.gz
	bawk '$$5=="$*" || $$5~/$*.intersec/ {print $$1,$$2,$$3,$$4";"$$5,$$5,"pos"}' $< > $@

%_neg.bed: ChIRP.bed.split.gz %_onlypos.bed
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

ALL_ssRNA=$(shell cat selected_ssRNA)

ALL_neg_pos_rand.bed: $(addsuffix .neg_pos_rand.bed,$(ALL_ssRNA))
	cut -f -4 $^ > $@
ALL_neg_pos.bed: $(addsuffix .neg_pos.bed,$(ALL_ssRNA))
	cut -f -4 $^ > $@

selected_ssRNA:
	@echo 'Human lncRNAs with idr peaks number > 100'
	@echo 'See this table: https://docs.google.com/spreadsheets/d/1iQGctC1ldu1oTwmaLEsrs8cB0ik3JWYCrBkvoDKZDpk/edit#gid=0'

pippo:
	echo $(GENCODE_DIR)

lncSmad7_onlypos.bed: /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/overlap.conservative.qTh005.mouse.regionPeak.gz
	bawk '$$5=="lncSmad7" || $$5~/lncSmad7.intersec/ {print $$1,$$2,$$3,$$4";"$$5,$$5,"pos"}' $< > $@
