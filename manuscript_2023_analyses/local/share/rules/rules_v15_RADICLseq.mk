##################
#  Generic rules

%.header_added: %
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@
%.header_added.gz: %.gz
	(bawk -M $< | cut -f 2 | transpose; zcat $< ) | gzip > $@
%.xlsx: %.gz
	zcat $< | tab2xlsx > $@
%.xlsx: %
	cat $< | tab2xlsx > $@


##############
#  AUC cmp

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.0_8_1_20_70_off_3.gz:
	matrix_reduce -t 'tpx_paramspace/*_ss0_unpairedWindow/*.neg_pos_rand.bed/min_length~8/max_length~-1/error_rate~20/guanine_rate~70/filter_repeat~off/consecutive_errors~3/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.gz' | \
	bawk 'NR==1 || $$2!="Duplex_ID" {split($$1,a,";"); print a[1],$$0}' | cut -f1,3- | gzip > $@
raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.10_10_1_20_10_off_3.gz:
	matrix_reduce -t 'tpx_paramspace/*_ss10_unpairedWindow/*.neg_pos_rand.bed/min_length~10/max_length~-1/error_rate~20/guanine_rate~10/filter_repeat~off/consecutive_errors~3/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.gz' | bawk 'NR==1 || $$2!="Duplex_ID" {split($$1,a,";"); print a[1],$$0}' | cut -f1,3- | gzip > $@
raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.%.AUC_cmp.tsv: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.%.gz
	zcat $< | ../../local/src/ROC.R neg_pos Stability_norm_undercount t_pot_norm Stability_best Score_best -d "<" -O raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.$*.roc > $@
