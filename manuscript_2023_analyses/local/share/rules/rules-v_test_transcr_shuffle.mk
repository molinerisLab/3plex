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

%_shuffle.fa: %.fa
	$(CONDA_ACTIVATE) /home/imoliner/.conda/envs/meme; \
	fasta-shuffle-letters -kmer $(FASTA_SHUFFLE_K) $< | sed 's/shuf/shuffle/' > $@
%.fa: $(VERSION_ssRNA_FASTA)/%.fa
	cp -a $< $@
%_posneg.bed: $(VERSION_BED)/%.neg_pos_rand.bed
	cp -a $< $@
%_posneg.fa: %_posneg.bed
	bedtools getfasta -fi $(GENOME_FA) -bed $< -name -fo $@
%_posneg.fasim.fa: %_posneg.bed
	bawk '{split($$4,a,";"); print $$1~3,a[2]"|"$$4,$$5,$$6}' $< | bedtools getfasta -name+ -fi $(GENOME_FA) -bed - | sed 's/::/|/' > $@

################
# TPX summary

#r_10_10_1_20_10_off_3
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
	cut -f4,6 $< | translate -a -v -e 0 <(bawk '{print $$1,$$14,$$15}' $^2) 1 | \
	bawk 'BEGIN{print "pos_neg","pred1","pred2"}{print $$4,$$2,$$3}' | gzip > $@
%.triplexAligner.summary.clean.gz: %_posneg.bed %.triplexAligner.summary.gz
	cut -f4,6 $< | translate -a -v -e 0 <(bawk 'NR>1{print $$12,$$7,$$10}' $^2 | find_best 1 3) 1 | \
	bawk 'BEGIN{print "pos_neg","pred1","pred2"} {print $$4,$$2,$$3}' | gzip > $@
%.fasimLongtarget.summary.clean.gz: %_posneg.bed %.fasimLongtarget.summary.gz
	cut -f4,6 $< | translate -a -v -e 0 <(bawk 'NR>1{print $$6,$$13,$$9}' $^2 | find_best 1 3) 1 | \
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

CDKN2B_AS1%_posneg.bed:
	cp $(VERSION_BED)/CDKN2B-AS1$*.neg_pos_rand.bed $@;
	sed -i 's/CDKN2B-AS1/CDKN2B_AS1/g' $@

CDKN2B_AS1.fa:
	cp $(VERSION_ssRNA_FASTA)/CDKN2B-AS1.fa $@;
	sed -i 's/CDKN2B-AS1/CDKN2B_AS1/g' $@