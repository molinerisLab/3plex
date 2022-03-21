ChIRP.bed.split.gz: ../../local/share/data/ReChIRP/idr_overlap_top1000/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.top1000.gz
	zcat $< | bedtools sort | gzip > $@

%.neg_pos.bed: %_onlypos.bed %_neg.bed
	cat $^ > $@

%_onlypos.bed: ChIRP.bed.split.gz
	bawk '$$5=="$*" || $$5~/$*.intersec/ {print $$1,$$2,$$3,$$4";"$$5,$$5,"pos"}' $< > $@

# idr + overlap top 1000 are considered positive but all the called peaks are considered in the blacklist
%_ALLpos.bed: ../../local/share/data/ReChIRP/idr_overlap_top1000/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.gz
	 bawk '$$1=="$*" {print $$2~4,"chirp_peak_" NR";"$$1,$$1,"pos"}' $< > $@

%_neg.bed: ChIRP.bed.split.gz %_onlypos.bed
	bawk '{print $$1,$$2,$$3,$$4";"$$5,"$*"}' $< > $@.tmp
	bedtools intersect -a $@.tmp -b $^2 -v \
	| bawk '{print $$0, "neg" }'> $@
	rm $@.tmp

%.neg_rand.excl.bed: $(BIOINFO_ROOT)/task/gencode/dataset/hsapiens/32/hg38.shuffle_blacklist.bed $(BIOINFO_ROOT)/task/gencode/dataset/hsapiens/32/gap.bed %_ALLpos.bed
	cut -f -3 $^ |  bedtools sort | bedtools merge > $@

%.neg_rand.bed: %_onlypos.bed $(BIOINFO_ROOT)/task/gencode/dataset/hsapiens/32/chrom.info.no_alt %.neg_rand.excl.bed
	bedtools shuffle -i $< -g $^2 -excl $^3 | bawk '{$$4="rand_"$$4; $$6="neg"; print}' > $@


%.neg_pos_rand.bed: %_onlypos.bed %.neg_rand.bed
	cat $^ > $@

ALL_ssRNA=$(shell cat selected_ssRNA_id)

ALL_neg_pos_rand.bed: $(addsuffix .neg_pos_rand.bed,$(ALL_ssRNA))
	cut -f -4 $^ > $@
ALL_neg_pos.bed: $(addsuffix .neg_pos.bed,$(ALL_ssRNA))
	cut -f -4 $^ > $@

selected_ssRNA:
	@echo 'Human lncRNAs with idr peaks number > 100'
	@echo 'See this table: https://docs.google.com/spreadsheets/d/1iQGctC1ldu1oTwmaLEsrs8cB0ik3JWYCrBkvoDKZDpk/edit#gid=0'
