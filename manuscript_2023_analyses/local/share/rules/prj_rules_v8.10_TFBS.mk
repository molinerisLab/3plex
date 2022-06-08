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

%.neg_rand.bed: %_onlypos.bed /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/37/chrom.info.no_alt %.neg_rand.excl.bed
	bedtools shuffle -i $< -g $^2 -excl $^3 | bawk '{$$4="rand_peak_"NR;$$5="neg"; print}' > $@

%.neg_pos_rand.bed: %_onlypos.bed %.neg_rand.bed
	cat $^ > $@

%.neg_pos_rand.bed.fimo_score: %_fimo/fimo.tsv %.neg_pos_rand.bed
	translate -a -f 3 -v -e 0 <(unhead $< | sort -k3,3 | unhead | grep -v '#' | find_best 3 7) 4 < $^2 > $@

.META: *.neg_pos_rand.bed.fimo_score
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


.SECONDARY: 

### TF binding sites ###
#
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1298-9#Sec23
# Supplementary file:
# https://static-content.springer.com/esm/art%3A10.1186%2Fs12859-016-1298-9/MediaObjects/12859_2016_1298_MOESM5_ESM.xls

ALL_TFs=CTCF_K562 E2F1_K562 GATA2_K562 IRF1_K562 MAX_K562 NFYA_K562 TAL1_K562 YY1_K562

E2F1_K562.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF326WMF/@@download/ENCFF326WMF.bed.gz'

IRF1_K562.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF019LOH/@@download/ENCFF019LOH.bed.gz'

TAL1_K562.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF665HHH/@@download/ENCFF665HHH.bed.gz'

CTCF_K562.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF041ODC/@@download/ENCFF041ODC.bed.gz'

MAX_K562.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF034DOZ/@@download/ENCFF034DOZ.bed.gz'

YY1_K562.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF193FQT/@@download/ENCFF193FQT.bed.gz'

GATA2_K562.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF930ZYR/@@download/ENCFF930ZYR.bed.gz'

NFYA_K562.idr_conservative.narrowPeaks.gz:
	wget -O $@ 'https://www.encodeproject.org/files/ENCFF195MIZ/@@download/ENCFF195MIZ.bed.gz'

IRF1_K562.meme:
	@echo 'Remove mouse motif!'
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0050/versions.meme'
GATA2_K562.meme:
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0036/versions.meme'
YY1_K562.meme:
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0095/versions.meme'
E2F1_K562.meme:
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0024/versions.meme'
CTCF_K562.meme:
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0139/versions.meme'
STAT1_K562.meme:
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0137/versions.meme'
MAX_K562.meme:
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0058/versions.meme'
NFYA_K562.meme:
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0060/versions.meme'
TAL1_K562.meme:
	@echo 'This motif is defined on Jaspar as TAL1::TCF3'
	wget -O $@ 'https://jaspar.genereg.net/api/v1/matrix/MA0091/versions.meme'

ALL_TFs.fimo_score_qval.AUC_matrix: $(addsuffix .neg_pos_rand.bed.fimo_score.AUC, $(ALL_TFs))
	matrix_reduce -t -l '$^' '*.neg_pos_rand.bed.fimo_score.AUC' | grep -v 'pred_1' > $@

ALL_TFs.fimo_score.AUC_tab:$(addsuffix .neg_pos_rand.bed.fimo_score.AUC, $(ALL_TFs))
	matrix_reduce -t -l '$^' '*.neg_pos_rand.bed.fimo_score.AUC' | grep -v 'pred_1' | bawk 'BEGIN{print "TF_cell_line","fimo_score_AUC"}{print $$1,$$5}' > $@

ALL_TFs.fimo_score.matrix.gz: $(addsuffix .neg_pos_rand.bed.fimo_score, $(ALL_TFs))
	matrix_reduce -t -l '$^' '*.neg_pos_rand.bed.fimo_score' | gzip > $@

.META: ALL_TFs.fimo_score.matrix.gz
	1	Gene_ID
	2	chr
	3	start
	4	stop
	5	peak_name
	6	motif_id
	7	motif_alt_id
	8	motif_start
	9	motif_stop
	10	strand
	11	score
	12	p_value
	13	q_value
	14	matched_sequence
	15	neg_pos

