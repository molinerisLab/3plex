##########################
#
#     General param
#

GENCODE_DIR=$(BIOINFO_REFERENCE_ROOT)/gencode/dataset/hsapiens/32
NCPU?=28


###########################
# 
#      Docker param
#

DOCKER_GROUP?=$(shell stat i-c %g $(PWD))
DOCKER_DATA_DIR?=$(PRJ_ROOT)



#######################################
#
#     Generic rules
#

%.header_added.bed: %.bed
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@
%.header_added: %
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@
%.header_added.gz: %.gz
	(bawk -M $< | cut -f 2 | transpose; zcat $< ) | gzip > $@
%.xlsx: %.gz
	zcat $< | tab2xlsx > $@
%.xlsx: %
	cat $< | tab2xlsx > $@
%.slop_tpx.bed: %.bed $(GENCODE_DIR)/chrom.info
	bedtools sort -i $< | bedtools slop -l $(TPX_DISTANCE) -r $(TPX_DISTANCE) -i stdin -g $^2 > $@
%.merged.bed: %.bed
	bedtools sort < $< | bedtools merge -c 4 -o distinct > $@
%.bb: %.bed $(GENCODE_DIR)/chrom.info
	bedToBigBed $< <(unhead $^2) $@
%.complement.bed: %.bed chrom.info.sorted
	bedtools complement -i <(bsort -k1,1V -k2,2n $<) -g $^2 > $@
chrom.info.sorted: $(GENCODE_DIR)/chrom.info
	unhead $< | bsort -k1,1V | grep -v _ > $@


######################################
#
#    LongTarget
#

cCRE.bed: $(BIOINFO_REFERENCE_ROOT)/encode-screen/dataset/v13/hg38-ccREs.bed
	ln -s $< $@

%.bed.fa: $(GENCODE_DIR)/GRCh38.primary_assembly.genome.fa %.bed
	bedtools getfasta -name -fi $< -bed $^2 -fo $@

%.mask.fa: $(GENCODE_DIR)/GRCh38.primary_assembly.genome.fa %.bed
	bedtools maskfasta -fi $< -fo $@ -bed $^2

chrs_mask/chr%: cCRE.complement.mask.fa
	mkdir -p `dirname $@`
	split_fasta -e -p $< 1 `dirname $@`
	for i in `dirname $@`/chr*; do perl -lpe '$$_.="|NA" if $$.==1' < $$i > $$i.tmp; mv $$i.tmp $$i; done

#all_lncRNA.LongTarget_mask.done: selected_ssRNA_id
#	for ssRNA in selected_ssRNA_id; do mkdir -p $$ssRNA; cd $$ssRNA; echo `seq 1 22` X Y M | tr " " "\n" | parallel -j 22 'LongTarget -f1 chr{} -f2 ssRNA.fa -ni 10 -ds 8 -lg 10'

%/ssRNA.fa: longest_transcripts.fa
	mkdir -p `dirname $@`
	get_fasta -i $* < $^ > $@


CHR?=chr1
#%/$(CHR)--$(CHR).t-TFOsorted: chrs_mask/$(CHR) chrs_mask/ssRNA.fa
#        cd chrs_mask;\
#        LongTarget -f1 $(CHR) -f2 ssRNA.fa -ni 10 -ds 8 -lg 10

%.ALL_LongTarget_mask.done: %/ssRNA.fa chrs_mask/chr1
	cd $*;\
	echo `seq 1 22` X Y M | tr " " "\n" | parallel -j 22 'LongTarget -f1 ../chrs_mask/chr{} -f2 ssRNA.fa -ni 10 -ds 8 -lg 10'
	touch $@

