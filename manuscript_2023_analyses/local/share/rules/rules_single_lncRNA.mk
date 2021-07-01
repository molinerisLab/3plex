lncRNA=$(shell basename $$PWD)

tpx: ../../hg38-ccREs.bed.fa.tpx ../../hg38-ccREs.bed
	bawk '$$sequence_id=="$(lncRNA)"' $< \
	| translate <(cut -f 4,5 $^2) 4 > $@

best_score.tpx: tpx
	find_best 4 7 < $< > $@

.META: *.tpx
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate

best_score.tpx.pos_neg: ../../hg38-NPC_H9.neg_pos.bed best_score.tpx
	cut -f 4- $< | translate -a -f 4 -v -e 0 $^2 1 > $@

best_score.tpx.pos_neg_no0: ../../hg38-NPC_H9.neg_pos.bed best_score.tpx
	translate -a -r -k <(cut -f 4- $<) 4 < $^2 > $@

.META: *.tpx.pos_neg best_score.tpx.pos_neg_no0
	1	sequence_id
	2	TFO_start
	3	TFO_end
	4	Duplex_ID
	5	TTS_start
	6	TTS_end
	7	Score
	8	Error_rate
	9	Errors
	10	Motif
	11	Strand
	12	Orientation
	13	Guanine_rate
	14	region_type
	15	pos_neg

