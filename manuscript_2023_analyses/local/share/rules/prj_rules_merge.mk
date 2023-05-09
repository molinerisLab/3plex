ChIRP.bed.split: ../../local/share/data/ChIRP_narrowPeaks.bed
	tr ";" "\t" < $< | grep -v CDKN2B-AS1 | bedtools sort > $@

ChIRP.bed.split.merge: ChIRP.bed.split
	cut -f 1-3,5 $< | union --allow-duplicates | bawk '{print $$1~3,"merged_ChIRP_"NR, $$4}'> $@

%_onlypos.bed: ChIRP.bed.split.merge
	bawk '$$5~/\<$*\>/ {print $$1,$$2,$$3,$$4,$$5,"pos"}' $< > $@

%.neg_pos.bed: %_onlypos.bed %_neg.bed
	cat $^ > $@

%_neg.bed: ChIRP.bed.split.merge %_onlypos.bed
	bawk '{print $$1,$$2,$$3,$$4,"$*"}' $< > $@.tmp
	bedtools intersect -a $@.tmp -b $^2 -v \
	| bawk '{print $$0, "neg" }'> $@
	rm $@.tmp

%.neg_rand.excl.bed: /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/hg38.shuffle_blacklist.bed %_onlypos.bed
	cut -f -3 $^ |  bedtools sort | bedtools merge > $@

%.neg_rand.bed: %_onlypos.bed /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/chrom.info.no_alt %.neg_rand.excl.bed
	bedtools shuffle -i $< -g $^2 -excl $^3 | bawk '{$$4="rand_"$$4; $$6="neg"; print}' > $@


%.neg_pos_rand.bed: %_onlypos.bed %.neg_rand.bed
	cat $^ > $@

ALL_ssRNA=$(shell cat selected_ssRNA_id)

ALL_neg_pos_rand.bed: $(addsuffix .neg_pos_rand.bed,$(ALL_ssRNA))
	cut -f -4 $^ > $@
ALL_neg_pos.bed: $(addsuffix .neg_pos.bed,$(ALL_ssRNA))
	cut -f -4 $^ > $@
