TRIPLEXATOR_PARAM?=-l 10 -e 20 -E 1
DOCKER_GROUP:=10001

##########################################3
#
#	triplexator
#


%.tfo: %.fa
	ssh epigen \
	sudo su - edoardo \
	/home/edoardo/src/triplexator/bin/triplexator $(TRIPLEXATOR_PARAM) -fm 1 -od . -o $@ -po -rm 2 -p $(CORES) -ss $<
%.tfo.log: %.tfo
	@echo done
%.tfo.summary: %.tfo
	@echo done
#%fa.tpx: ssRNA.fa %fa
#	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) triplexator:v0.02 bash -c "cd $(PWD); triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 0     -o $@ -rm 2 -p $(CORES) -ss $< -ds $^2" 
#see rules.mk
%fa.tpx_aln: ssRNA.fa %fa
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) triplexator:v0.02 bash -c "cd $(PWD); triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 1 -po -o $@ -rm 2 -p $(CORES) -ss $< -ds $^2" 

.META: *fa.tpx
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

%.fa.tpx.around_validated_tfo: %.fa.tpx
	bawk '{named_tfo="no_tpx"; \
		if($$Duplex_ID!="no_tpx"){ \
			named_tfo="no_name";\
			if ($$TFO_start>=83 && $$TFO_start<=89) {named_tfo="86"} \
			if($$TFO_start>=45 && $$TFO_start<=51)  {named_tfo="48"} \
		} \
		print $$0,named_tfo} \
	' $< > $@

.META: *fa.around_validated_tfo
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
	14	tfo_name

#########################################################
#
#	from tpx to bed (genome coords)
#

%.fa.tpx.tts_genom_coords.pre: %.fa.tpx %
	unhead $< | translate -r -a -f 4 <(cut -f -4 $^2) 4 > $@
%.fa.tpx.around_validated_tfo.tts_genom_coords.pre: %.fa.tpx.around_validated_tfo %
	unhead $< | translate -r -a -f 4 $^2 4 > $@

.META: *.tpx.tts_genom_coords.pre
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
	14	region_chr
	15	region_b
	16	region_e

.META: *.tpx.tts_genom_coords
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
	14	region_chr
	15	region_b
	16	region_e
	17	TTS_start_genom
	18	TTS_end_genom

.META: Srsf3_region.bed.fa.tpx.around_validated_tfo.tts_genom_coords
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
	14	tfo_name
	15	region_chr
	16	region_b
	17	region_e
	18	TTS_start_genom
	19	TTS_end_genom

.META: *.tpx.around_validated_tfo.tts_genom_coords.pre *.tpx.around_validated_tfo.tts_genom_coords
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
	14	tfo_name
	15	region_chr
	16	region_b
	17	region_e
	18	TTS_start_genom
	19	TTS_end_genom

%.tpx.tts_genom_coords: %.tpx.tts_genom_coords.pre
	bawk '{TTS_start_genom=$$TTS_start+$$region_b; TTS_end_genom=$$TTS_end+$$region_b; print $$0,TTS_start_genom,TTS_end_genom}' $< > $@
%.tpx.around_validated_tfo.tts_genom_coords: %.tpx.around_validated_tfo.tts_genom_coords.pre
	bawk '{TTS_start_genom=$$TTS_start+$$region_b; TTS_end_genom=$$TTS_end+$$region_b; print $$0,TTS_start_genom,TTS_end_genom}' $< > $@

%.tpx.tts_genom_coords.bed: %.tpx.tts_genom_coords
	bawk '{print $$region_chr,$$TTS_start_genom,$$TTS_end_genom, $$sequence_id "-" $$Duplex_ID, $$Score,$$Strand,$$sequence_id,$$TFO_start,$$TFO_end,$$Error_rate,$$Errors,$$Motif,$$Orientation,$$Guanine_rate,$$Duplex_ID,$$region_b,$$region_e}' $< \
	| id2count -a 4  >$@
%.tpx.tts_genom_coords.aln.bed: %.tpx.tts_genom_coords.aln
	bawk '{print $$region_chr,$$TTS_start_genom,$$TTS_end_genom, $$sequence_id "-" $$Duplex_ID, $$Score,$$Strand,$$sequence_id,$$TFO_start,$$TFO_end,$$Error_rate,$$Errors,$$Motif,$$Orientation,$$Guanine_rate,$$Duplex_ID,$$region_b,$$region_e,$$alignment1,$$alignment2,$$alignment3,$$alignment4}' $< \
	| id2count -a 4  >$@
input_regions.bed.fa.tpx.tts_genom_coords.minimal.bed: input_regions.bed.fa.tpx.tts_genom_coords.bed
	cut -f -5 $< | bedtools sort | uniq > $@
%.tpx.around_validated_tfo.tts_genom_coords.bed: %.tpx.around_validated_tfo.tts_genom_coords
	bawk '{print $$region_chr,$$TTS_start_genom,$$TTS_end_genom, $$sequence_id "-" $$tfo_name "-" $$Duplex_ID, $$Score,$$Strand,$$sequence_id,$$TFO_start,$$TFO_end,$$Error_rate,$$Errors,$$Motif,$$Orientation,$$Guanine_rate,$$tfo_name,$$Duplex_ID,$$region_b,$$region_e}' $< \
	| id2count -a 4 >$@

.META: *.tpx.around_validated_tfo.tts_genom_coords.bed
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	tfo_name
	16	region_name
	17	region_b
	18	region_e

.META: *.tpx.tts_genom_coords.bed
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	region_name
	16	region_b
	17	region_e


.META: *.tts_genom_coords.regulated_genes.bed
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	region_name
	16	region_b
	17	region_e
	18	gene_chr
	19	TSS_or_PIR_b
	20	TSS_or_PIR_e
	21	GeneID
	22	TSS_or_PIR

##########################################################################
#
#	repeat in tts
#

tts.simple_repeat: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed /sto1/ref/bioinfotree/task/ucsc-tracks/dataset/mm10/repeat_simple.extended.bed.gz
	cut -f -4 $< | bedtools intersect -a - -b $^2 -loj | cut -f 4,8 > $@

cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.repeat_simple: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed /sto1/ref/bioinfotree/task/ucsc-tracks/dataset/mm10/repeat_rmsk.simple.bed.gz
	bedtools intersect -a $< -b $^2 -loj > $@
cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.repeat_noSimple: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed /sto1/ref/bioinfotree/task/ucsc-tracks/dataset/mm10/repeat_rmsk.noSimple.bed.gz
	bedtools intersect -a $< -b $^2 -loj > $@

.META: *.tpx.tts_genom_coords.bed.repeat_noSimple *.tpx.tts_genom_coords.bed.repeat_simple
	1	region_chr
	2	TTS_start
	3	TTS_end
	4	TTS_name
	5	Score
	6	Strand
	7	sequence_id
	8	TFO_start
	9	TFO_end
	10	Error_rate
	11	Errors
	12	Motif
	13	Orientation
	14	Guanine_rate
	15	region_name
	16	region_b
	17	region_e
	18	repeat_chr
	19	repeat_b
	20	repeat_e
	21	repeat_name

tts_simple_repeat_enric: cCRE.slop_tpx.merged.bed.fa.tpx.tts_genom_coords.bed.repeat_simple
	parallel 'bawk '\''$$Score>{} {print $$1~3,$$repeat_name}'\'' $< | bedtools sort | bedtools merge -c 4 -o distinct | cut -f 4 | symbol_count -n | bsort -k2,2nr | bawk '\''$$2>=0.001'\'' | append_each_row -B {}' ::: 13 16 19 21 | tab2matrix -r repeat -t > $@

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
hg38-ccREs.bed.fa.tpx.%_pos_neg.fischer_select_cutoff: hg38-ccREs.bed.fa.tpx.%_pos_neg
	bawk '{print $$1";"$$4,$$3,$$5}' $< | ./fischer_select_cutoff.py -a greater | tr ";" "\t" > $@

.META: hg38-ccREs.bed.fa.tpx.*_pos_neg.fischer_select_cutoff
	1	lncRNA
	2	cCRE
	3	score
	4	greater_positive
	5	greater_negative
	6	lower_positive
	7	lower_negative
	8	oddsratio
	9	pvalue

hg38-ccREs.bed.fa.tpx.ALL_pos_neg: hg38-ccREs.bed.fa.tpx.NPC_pos_neg hg38-ccREs.bed.fa.tpx.H9_pos_neg
	matrix_reduce 'hg38-ccREs.bed.fa.tpx.*_pos_neg' -l '$^' | fasta2tab > $@

.META: hg38-ccREs.bed.fa.tpx.ALL_pos_neg
	1	cCRE_source	NPC
	2	lncRNA		LINC00461
	3	cCRE_id		EH38E0065969
	4	score		10
	5	cCRE_type	dELS
	6	pos_neg		neg

hg38-ccREs.bed.fa.tpx.ALL_pos_neg.fischer_select_cutoff: hg38-ccREs.bed.fa.tpx.NPC_pos_neg.fischer_select_cutoff hg38-ccREs.bed.fa.tpx.H9_pos_neg.fischer_select_cutoff
	matrix_reduce 'hg38-ccREs.bed.fa.tpx.*_pos_neg.fischer_select_cutoff' -l '$^' | fasta2tab | bawk '{print $$2~10,$$1}' > $@

.META: *.fischer_select_cutoff
	1	lncRNA
	2	cCRE
	3	score
	4	greater_positive
	5	greater_negative
	6	lower_positive
	7	lower_negative
	8	oddsratio
	9	pvalue


hg38-cCRE_specific_H9-NPC: hg38-H9.neg_pos.bed hg38-NPC.neg_pos.bed
	matrix_reduce 'hg38-*.neg_pos.bed' -l '$^' | fasta2tab | bawk '$$7=="pos" {print $$5,$$6,$$1}' | sed 's/NPC/NPC/' | collapsesets 3 | grep -v ';' > $@
hg38-ccREs.bed.fa.tpx.best_score.specific_H9-NPC: hg38-cCRE_specific_H9-NPC hg38-ccREs.bed hg38-ccREs.bed.fa.tpx.best_score
	translate -a -r -k -d <(cat $< | translate -f 2 <(cut -f 4,5 $^2) 1) 2 < $^3 \
	| grep -v ';' > $@                      *perdiamo i ccRE che cambiano tipo ad H9 a NPC, ma sono pochi e nessun dELS*
hg38-ccREs.bed.fa.tpx.best_score.specific_H9-NPC.fischer_select_cutoff: hg38-ccREs.bed.fa.tpx.best_score.specific_H9-NPC
	bawk '{if($$5=="NPC"){$$5="pos"}else{$$5="neg"} print $$1";"$$4,$$3,$$5}' $< | fischer_select_cutoff  | tr ";" "\t" > $@
