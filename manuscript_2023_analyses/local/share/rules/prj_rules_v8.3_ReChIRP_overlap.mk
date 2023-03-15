3plex/TERC/TERC.neg_pos_rand.bed/tpx.all_scores.gz: ../best_single_params.triplex_ssRNA_scores.header_added.gz
	mkdir -p `$@`; \
	bawk 'NR==1{print} NR>1 && $1=="TERC"{print}' $< | gzip > $@
3plex/TERC/TERC.neg_pos_rand.bed/tpx.specific_score.gz: 3plex/TERC/TERC.neg_pos_rand.bed/tpx.all_scores.gz
	zcat $< | bawk '{print $$2,$$18,$$4}' | gzip > $^2
