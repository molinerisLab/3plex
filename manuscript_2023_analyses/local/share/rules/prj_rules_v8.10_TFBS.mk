CTCF_GM12878.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF796WRU/@@download/ENCFF796WRU.bed.gz'

CTCF_GM12878.meme: 
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0139.1.meme'

STAT3_GM12878.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF719YDJ/@@download/ENCFF719YDJ.bed.gz'

STAT3_GM12878.meme:
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0144.2.meme'

%.fa: $(BIOINFO_ROOT)/task/gencode/dataset/hsapiens/37/GRCh38.primary_assembly.genome.clean_id.fa %.neg_pos_rand.bed
	bedtools getfasta -name -fi $< -bed $^2 > $@

%_fimo/fimo.tsv: %.fa %.meme
	conda activate meme; \
	fimo --oc `dirname $@` $^2 $<

%_onlypos.bed: %.idr_conservative.narrowPeaks.gz
	bawk '{print $$1,$$2,$$3,"chip_peak_"NR,"pos"}' $< > $@

%.neg_rand.excl.bed: $(BIOINFO_ROOT)/task/gencode/dataset/hsapiens/37/hg38.shuffle_blacklist.bed $(BIOINFO_ROOT)/task/gencode/dataset/hsapiens/37/gap.bed %_onlypos.bed
	cut -f -3 $^ |  bedtools sort | bedtools merge > $@

%.neg_rand.bed: %_onlypos.bed /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/chrom.info.no_alt %.neg_rand.excl.bed
	bedtools shuffle -i $< -g $^2 -excl $^3 | bawk '{$$4="rand_peak_"NR;$$5="neg"; print}' > $@

%.neg_pos_rand.bed: %_onlypos.bed %.neg_rand.bed
	cat $^ > $@

%.neg_pos_rand.bed.fimo_score: %_fimo/fimo.tsv %.neg_pos_rand.bed
	translate -a -f 3 -v -e 0 <(unhead $< | sort -k3,3 | unhead | grep -v '#' | find_best 3 7) 4 < $^2 > $@

.META: %.neg_pos_rand.bed.fimo_score
	1	chr
	2	start
	3	stop
	4	peak_name
	5	motif_id
	6	motif_alt_id
	7	motif_start
	8	motif_stop
	9	strand
	10	score
	11	p_value
	12	q_value
	13	matched_sequence
	14	neg_pos

%.neg_pos_rand.bed.fimo_score.AUC: %.neg_pos_rand.bed.fimo_score.header_added
	cat $< | ../../local/src/ROC.R neg_pos score q_value -d "<" -O '$*.auc_comparison_plot' > $@
