GENCODE_DIR=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/$(SPECIES)/$(VERSION)
CHROM_INFO=$(GENCODE_DIR)/chrom.info
GENOME_FA=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/hsapiens/32/GRCh38.primary_assembly.genome.clean_id.fa 
SEED?=42

CONDA_ROOT="/opt/conda"
CONDA_VERSION="miniconda3"
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate

THREADS=16

SAMPLES=$(shell cat selected_ssRNA)

.SECONDARY:

##################
#  Prepare Inputs

#%.fa: $(GENCODE_DIR)/transcripts.fa.gz
#	zcat $< | fasta2oneline | tr "|" "\t" | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,10 | tab2fasta | fold | get_fasta -i $* > $@
ChIRP.bed.split.gz: $(PEAKS)
	zcat $< | bedtools sort | tr "-" "_" | gzip > $@
%_pos.bed: ChIRP.bed.split.gz
	bawk '$$5=="$*" {print $$1,$$2,$$3,$$4";"$$5,"pos"}' $< > $@
%_neg.bed: $(addsuffix _pos.bed, $(SAMPLES))
	cat $^ | tr ";" "\t" | bawk '$$5!="$*" {$$4="rand_"$$4";$*"; $$5="neg"; print}' | cut -f-5 > $@
%_posneg.bed: %_pos.bed %_neg.bed
	bawk '{$$4=$$4";"$$5; print}' $< $^2 > $@
%_posneg.fa: %_posneg.bed
	bedtools getfasta -fi $(GENOME_FA) -bed $< -name -fo $@
%_posneg.fasim.fa: %_posneg.bed
	bawk '{split($$4,a,";"); print $$1~3,a[2]"|"$$4,$$5,$$6}' $< | bedtools getfasta -name+ -fi $(GENOME_FA) -bed - | sed 's/::/|/' > $@

################
# TPX summary

%.3plex.summary.gz: %.fa %_posneg.fa
	docker run -u `id -u`:`id -g` --rm -v $$PWD:$$PWD imolineris/3plex:v0.1.2-beta -j $(THREADS) -l 8 -L 1 -e 20 -s 0 -g 70 -c 3 $$PWD/$< $$PWD/$^2 $$PWD
	mv $*_ssmasked-$*_posneg.tpx.summary.gz $@
%.triplexAligner.summary.gz: %.fa %_posneg.fa
	docker run -u `id -u`:`id -g` --rm -v $$PWD:$$PWD triplex_aligner $$PWD/$^2 $$PWD/$< hs | gzip > $@
%.fasimLongtarget.summary.gz: %.fa %_posneg.fasim.fa
	docker run --rm -v $$PWD:$$PWD fasim -f1 $$PWD/$^2 -f2 $$PWD/$< -O $$PWD
	cat $*-$*-fastSim-TFOsorted | gzip > $@
	rm $*-$*-fastSim-TFOsorted

#########################
# Single summary clean 

%.3plex.summary.clean.gz: %_posneg.bed %.3plex.summary.gz
	cut -f4,5 $< | translate -a -v -e 0 <(bawk '{print $$1,$$14,$$15}' $^2) 1 | \
	bawk 'BEGIN{print "pos_neg","pred1","pred2"}{print $$4,$$2,$$3}' | gzip > $@
%.triplexAligner.summary.clean.gz: %_posneg.bed %.triplexAligner.summary.gz
	cut -f4,5 $< | translate -a -v -e 0 <(bawk 'NR>1{print $$12,$$7,$$10}' $^2 | find_best 1 3) 1 | \
	bawk 'BEGIN{print "pos_neg","pred1","pred2"} {print $$4,$$2,$$3}' | gzip > $@
%.fasimLongtarget.summary.clean.gz: %_posneg.bed %.fasimLongtarget.summary.gz
	cut -f4,5 $< | translate -a -v -e 0 <(bawk 'NR>1{print $$6,$$9,$$13}' $^2 | find_best 1 2) 1 | \
	bawk 'BEGIN{print "pos_neg","pred1","pred2"} {print $$4,$$2,$$3}' | gzip > $@

#####################
# All summary clean

3plex.summary.clean.gz: $(addsuffix .3plex.summary.clean.gz, $(SAMPLES))
	zcat $^ | bawk 'NR==1 || $$1!="pos_neg"' | gzip > $@
triplexAligner.summary.clean.gz: $(addsuffix .triplexAligner.summary.clean.gz, $(SAMPLES))
	zcat $^ | bawk 'NR==1 || $$1!="pos_neg"' | gzip > $@ 
fasimLongtarget.summary.clean.gz: $(addsuffix .fasimLongtarget.summary.clean.gz, $(SAMPLES))
	zcat $^ | bawk 'NR==1 || $$1!="pos_neg"' | gzip > $@ 

##########
# AUC cmp

%.summary.clean.AUC_cmp.tsv: %.summary.clean.gz
	$(CONDA_ACTIVATE) /home/cciccone/.conda/envs/pROC_Env; \
	zcat $< | ../../local/src/ROC.R pos_neg pred1 pred2 -d "<" -O $*.roc > $@
