CHROM_INFO=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/hsapiens/32/chrom.info
GENOME_FA=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/hsapiens/32/GRCh38.primary_assembly.genome.clean_id.fa 
SEED?=42

CONDA_ROOT="/opt/conda"
CONDA_VERSION="miniconda3"
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate

THREADS=16

SAMPLES=$(shell cat selected_ssRNA)

ChIRP.bed.split.gz: ../../local/share/data/ReChIRP/idr_overlap_top1000/all_reproducibility.idr_conservative-idr_optimal_peak-overalp_conservative-overlap_optimal.regionPeak.top1000.gz
	zcat $< | bedtools sort | gzip > $@
%_pos.bed: ChIRP.bed.split.gz
	bawk '$$5=="$*" {print $$1,$$2,$$3,$$4";"$$5,"pos"}' $< > $@

rand.excl.bed: /home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/32/hg38.shuffle_blacklist.bed /home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/32/gap.bed
	cut -f -3 $< $^2 | bedtools sort | bedtools merge > $@

%_neg.bed: rand.excl.bed %_pos.bed /home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/32/chrom.info.no_alt
	bedtools shuffle -excl $< -i $^2 -g $^3 -seed $(SEED) \
	| bawk '{$$4="rand_"$$4; $$5="neg"; print}' > $@

%_posneg.bed: %_pos.bed %_neg.bed
	bawk '{$$4=$$4";"$$5; print}' $< $^2 > $@

%.fa: ../v8.6_ReChIRP_idr_overlap_top1000/%.fa
	cp -a $< $@

%_posneg.fa: %_posneg.bed
	bedtools getfasta -fi $(GENOME_FA) -bed $< -name+ -fo $@

%.3plex.summary.gz: %.fa %_posneg.fa
	docker run -u `id -u`:`id -g` -it --rm -v $$PWD:$$PWD imolineris/3plex:v0.1.2-beta -j $(THREADS) -l 8 -L 1 -e 20 -s 0 -g 70 -c 3 $$PWD/$< $$PWD/$^2 $$PWD
	mv $*_ssmasked-$*_posneg.tpx.summary.gz $@

%.triplexAligner.summary.gz: %.fa %_posneg.fa
	docker run -u `id -u`:`id -g` -it --rm -v $$PWD:$$PWD triplex_aligner $$PWD/$^2 $$PWD/$< nonesiste | gzip > $@

%.fasimLongtarget.summary.gz: %.fa %_posneg.fasim.fa
	docker run --rm -v $$PWD:$$PWD fasim -f1 $$PWD/$^2 -f2 $$PWD/$< -O $$PWD
	cat $*-$*-fastSim-TFOsorted | gzip > $@
	rm $*-$*-fastSim-TFOsorted

%.3plex.summary.clean.gz: %.3plex.summary.gz
	bawk 'NR==1 {$$1="Duplex_ID;lncRNA;pos_neg"; $$14="pred1"; $$15="pred2"} {print}' $< | tr ";" "\t" | cut -f 3,16,17 | gzip > $@

3plex.summary.clean.gz: $(addsuffix .3plex.summary.clean.gz, $(SAMPLES))
	zcat $^ | bawk 'NR==1 || $$1!="pos_neg"' | gzip > $@

triplexAligner.summary.clean.gz: $(addsuffix .3plex.summary.clean.gz, $(SAMPLES))
	# todo

%.fasimLongtarget.summary.clean.gz: %.fasimLongtarget.summary.gz
	bawk 'NR==1 {$$6="Duplex_ID;lncRNA;pos_neg"; $$9="pred1"; $$13="pred2"} {print $$6,$$9,$$13}' $< | tr ";" "\t" | cut -f3- | gzip > $@

%.summary.clean.AUC_cmp.tsv: %.summary.clean.gz
	$(CONDA_ACTIVATE) /home/cciccone/.conda/envs/pROC_Env; \
	zcat $< | ../../local/src/ROC.R pos_neg pred1 pred2 > $@

#%_posneg.fasim.fa: %_posneg.bed
#	bedtools getfasta -fi $(GENOME_FA) -bed $< -name+ | sed 's/::/|/; s/:/|/;' > $@
%_posneg.fasim.fa: %_posneg.bed
	bawk '{split($$4,a,";"); print $$1~3,a[2]"|"$$4,$$5,$$6}' $< | bedtools getfasta -name+ -fi $(GENOME_FA) -bed - | sed 's/::/|/' > $@
