GENCODE_DIR=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/$(SPECIES)/$(VERSION)

SEED?=42

CONDA_ROOT="/opt/conda"
CONDA_VERSION="miniconda3"
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate

THREADS=16

SAMPLES=$(shell cat selected_ssRNA)

.SECONDARY:

##################
#  Prepare Inputs

%.fa: $(VERSION_ssRNA_FASTA)/%.fa
	cp -a $< $@
ChIRP.bed.split.gz: $(PEAKS)
	zcat $< | bedtools sort | tr "-" "_" | gzip > $@
%_pos.bed: ChIRP.bed.split.gz
	bawk '$$5=="$*" {print $$1,$$2,$$3,$$4";"$$5,"pos"}' $< > $@
rand.excl.bed: $(GENCODE_DIR)/$(GENOME).shuffle_blacklist.bed $(GENCODE_DIR)/gap.bed
	cut -f -3 $< $^2 | bedtools sort | bedtools merge > $@
%_neg.bed: rand.excl.bed %_pos.bed $(GENCODE_DIR)/chrom.info.no_alt
	bedtools shuffle -excl $< -i $^2 -g $^3 -seed $(SEED) \
	| bawk '{$$4="rand_"$$4; $$5="neg"; print}' > $@
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
	cut -f4,5 $< | translate -a -v -e 0 <(bawk 'NR>1{print $$6,$$13,$$9}' $^2 | find_best 1 3) 1 | \
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

AUC_singles.tsv: $(addsuffix .3plex.summary.clean.AUC_cmp.tsv, $(SAMPLES)) $(addsuffix .triplexAligner.summary.clean.AUC_cmp.tsv, $(SAMPLES)) $(addsuffix .fasimLongtarget.summary.clean.AUC_cmp.tsv, $(SAMPLES))
	matrix_reduce -t -l '$^' '*\.*\.summary.clean.AUC_cmp.tsv' | tr ";" "\t" | cut -f1,2,4,6 | grep -v AUC > $@
AUC_singles.matrix.xlsx: AUC_singles.tsv
	cat $< | cut -f1,2,4 | tab2matrix -r ssRNA | tab2xlsx > $^2

NEAT1.triplexAligner.summary.gz.time:
	time remake NEAT1.triplexAligner.summary.gz > $@
