REGION_TYPE=$(shell basename $$PWD)

tpx.pos_neg: ../../best_score.tpx.pos_neg
	bawk '$$region_type=="$(REGION_TYPE)"' $< > $@
tpx.pos_neg_no0: ../../best_score.tpx.pos_neg_no0
	bawk '$$region_type=="$(REGION_TYPE)"' $< > $@

.META: tpx.pos_neg tpx.pos_neg_no0
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

%.header_added: %
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@

tpx.%.bmodel: tpx.%.header_added
	bmodel -s -c -u -B $< "pos_neg~Score+Guanine_rate" > $@

tpx.%.bmodel_bootstrap: tpx.%.header_added
	bmodel -s -c -u -b 100 -B $< "pos_neg~Score+Guanine_rate" > $@ 2>$@.log
