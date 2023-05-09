ChIRP.bed.split.gz: ../../local/share/data/ReChIRP/mmusculus/idr_overlap_top1000/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.top1000.gz
	zcat $< | bedtools sort | gzip > $@

%.neg_pos.bed: %_onlypos.bed %_neg.bed
	cat $^ > $@

%_onlypos.bed: ChIRP.bed.split.gz
	bawk '$$5=="$*" || $$5~/$*.intersec/ {print $$1,$$2,$$3,$$4";"$$5,$$5,"pos"}' $< > $@

# idr + overlap top 1000 are considered positive but all the called peaks are considered in the blacklist
%_ALLpos.bed: ../../local/share/data/ReChIRP/mmusculus/idr_overlap_top1000/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.gz
	bawk '$$1=="$*" {print $$2~4,"chirp_peak_" NR";"$$1,$$1,"pos"}' $< > $@

%_neg.bed: ChIRP.bed.split.gz %_onlypos.bed
	bawk '{print $$1,$$2,$$3,$$4";"$$5,"$*"}' $< > $@.tmp
	bedtools intersect -a $@.tmp -b $^2 -v \
	| bawk '{print $$0, "neg" }'> $@
	rm $@.tmp

%.neg_rand.excl.bed: $(GENCODE_DIR)/mm10.shuffle_blacklist.bed $(GENCODE_DIR)/gap.bed %_ALLpos.bed
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
	@echo 'Mouse lncRNA with idr cons + idr opt + ovelap cons top 1000'
	@echo 'See this table: https://docs.google.com/spreadsheets/d/1iQGctC1ldu1oTwmaLEsrs8cB0ik3JWYCrBkvoDKZDpk/edit#gid=0'

pippo:
	echo $(GENCODE_DIR)

#lncSmad7_onlypos.bed: /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.top1000.gz
#	bawk '$$2=="lncSmad7" {print $$4,$$5,$$6,"chirp_peak_"NR";"$$2,$$2,"pos"}' $< > $@
3plex/lncSmad7/lncSmad7.neg_pos_rand.bed/tpx.all_scores.gz: ../best_single_params.triplex_ssRNA_scores.header_added.gz
	bawk 'NR==1{print} NR>1 && $$1=="lncSmad7"{print}' $< | gzip > $@
3plex/lncSmad7/lncSmad7.neg_pos_rand.bed/tpx.specific_score.gz: 3plex/lncSmad7/lncSmad7.neg_pos_rand.bed/tpx.all_scores.gz
	zcat $< | bawk '{print $$2,$$13,$$4}' | gzip > $@


tpx_paramspace_AUC_cmp.triplex_ssRNA.mean_AUC.gz: tpx_paramspace_AUC_cmp.gz selected_ssRNA
	zgrep -v TTS_covered_frac $< | bawk 'BEGIN{FS = "\t";OFS = ";"}{print $$1"\t"$$2,$$4~10"\t"$$12}' | sort -u | filter_1col 1 $^2 | cut -f2,3 | sort -k1,1 | stat_base -g -a | sort -k2,2nr | gzip > $^3

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.gz:
	matrix_reduce -t 'tpx_paramspace/*_ss*_unpairedWindow/*.neg_pos_rand.bed/min_length~*/max_length~*/error_rate~*/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.gz' | \
	grep -v 'Duplex_ID' | gzip > $@
