#layout
#cat ../../local/share/data/selected_lncRNA | parallel -j 1 --lb 'cd single_lncRNA; dataset_clone_recursive MIAT {} | sh'

###########################
#
#      gencode param
#

GENCODE_SPECIES?=hsapiens
GENCODE_VERSION?=32
GENCODE_GENOME?=GRCh38.primary_assembly
GENCODE_DIR=$(BIOINFO_ROOT)/task/gencode/dataset/$(GENCODE_SPECIES)/$(GENCODE_VERSION)

GENCODE_ANNOTATION?=basic.annotation
COUNT_REF_GTF?=$(GENCODE_DIR)/$(GENCODE_ANNOTATION).gtf
COUNT_REF_GTF_COMP?=$(GENCODE_DIR)/primary_assembly.annotation.gtf
COUNT_REF_BED=$(GENCODE_DIR)/$(GENCODE_ANNOTATION).bed

###########################
#
#   triplex param
#

TRIPLEXATOR_PARAM?=-l 10 -L -1 -e 20 -E -1
TDF_PARAM?=I -ccf 1 -l 10 -e 20 -par L_-1_-E_-1

cCRE?=$(BIOINFO_REFERENCE_ROOT)/encode-screen/dataset/v13/hg38-NPC.bed
#cCRE?=$(BIOINFO_REFERENCE_ROOT)/encode-screen/dataset/v13/GRCh38-ccREs.bed


# file fasta unico con tutti i lncs selezionati (da usare come input per triplexator):
ssRNA.fa: $(BIOINFO_REFERENCE_ROOT)/gencode/dataset/$(GENCODE_SPECIES)/$(GENCODE_VERSION)/gencode.v32.transcripts.fa.gz ../../local/share/data/selected_lncRNA
	zcat $< | fasta2oneline | tr "|" "\t" | filter_1col 6 <(cut -f 1 $^2) | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,10 | tab2fasta | fold > $@

# ciclo che itera in lista di gene id dei lncRNA selezionati per creare un file fasta separato per ciascuno (da usare come inpute per TDF)
#single_fasta_done: $(BIOINFO_REFERENCE_ROOT)/gencode/dataset/$(GENCODE_SPECIES)/$(GENCODE_VERSION)/gencode.v32.transcripts.fa.gz selected_lncRNA
#	while read i; do zcat $< | fasta2oneline | tr "|" "\t" | grep $$i  | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,10 | tab2fasta | fold > TDF/$$i.fa; done < $<; touch $@


%.fa: ssRNA.fa
	get_fasta -i $* < $< > $@

############################
#
#       region definition
#

cCRE.bed: $(cCRE)
	bawk '$$10=="PLS" || $$10=="dELS" || $$10=="pELS"' $< > $@
%.bed.fa: $(GENCODE_DIR)/GRCh38.primary_assembly.genome.clean_id.fa %.bed
	bedtools getfasta -name -fi $< -bed $^2 -fo $@

#peak.fa: /sto1/epigen/Fatemeh_hESC/dataset/PARCLIP_P300_v1/Fatemeh-v2-P300_merged_reps.gencode_exon_annotation.seq.header_aded.gz ../../local/share/data/selected_lncRNA
#	bawk 'NR>1 {print $$exon_geneID, $$peak_seq}' $< | filter_1col 1 $^2 | id2count -b 1 | bawk '{print ">"$$1; print $$2}' > $@

###########################
#       triplexator       #
###########################

%fa.tfo: %fa
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) triplexator:v0.02 bash -c "cd $(PWD); triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 0     -o $@ -rm 2 -p $(CORES) -ss $<"

%fa.tts: %fa
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) triplexator:v0.02 bash -c "cd $(PWD); triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 0     -o $@ -rm 2 -p $(CORES) -ds $<" 

%fa.tpx: ssRNA.fa %fa
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $(DOCKER_DATA_DIR):$(DOCKER_DATA_DIR) -v $(SCRATCH):$(SCRATCH) triplexator:v0.02 bash -c "cd $(PWD); triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 0     -o $@ -rm 2 -p $(CORES) -ss $< -ds $^2" 
	bawk '$$Score>=$(TRIPLEXATOR_SCORE_CUTOFF)' $@ >$@.tmp
	mv $@.tmp $@
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

%.xlsx: %.gz
	zcat $< | tab2xlsx > $@
%.xlsx: %
	cat $< | tab2xlsx > $@

.META: *.tfo
	1	SequenceID
	2	Start
	3	End
	4	Score
	5	Motif
	6	ErrorRate
	7	Errors
	8	GuanineRate
	9	Duplicates
	10	TFO
	11	Duplicate locations

ssRNA.fa.tfo.coverage: ssRNA.fa.tfo
	bawk '$$ErrorRate<=0.05 && $$Score>=13 {print $$SequenceID,$$Start,$$End}' $< | bedtools sort | bedtools merge | bawk '{print $$1,$$3-$$2}' | stat_base -g -t > $@
ssRNA.len: /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/gencode.v32.transcripts.fa.gz ../../local/share/data/selected_lncRNA
	zcat $< | fasta2oneline | tr "|" "\t" | filter_1col 6 <(cut -f 1 $^2) | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,7 > $@
ssRNA.fa.tfo.coverage.norm: ssRNA.len ssRNA.fa.tfo.coverage
	translate -a -r $< 1 < $^2 | bawk '{print $$0,$$2/$$3}' > $@

###################
#  	TDF       #
###################

CORES=10
RGTDATA_DIR=/sto1/ref/bioinfotree/task/rgtdata/0.13.1
DOCKER_GROUP=1000

TSS?=/sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/primary_assembly.annotation.tss.bed
CHROM_SIZE?=/sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/chrom.info
REGIONS?=GeneIDs_NPC_TD-up.tss2kb.bed


%/index.html: GeneIDs_NPC_TD-up.txt TDF/%.fa
	docker run -u `id -u`:$(DOCKER_GROUP) --rm -v $$PWD:$$PWD -v $(RGTDATA_DIR):/opt/rgtdata -e RGTDATA=/opt/rgtdata rgt:0.13.1 bash -c "cd $$PWD; rgt-TDF promotertest -r $^2 -de $< -rn $* -organism hg38 -o TDF/promoter_test/$*_l10e20E1 $(TDF_PARAM)" 

NPC-edger.toptable.gz: ../../../RNAseq/dataset/v1/DGE_H9_NPC/edger.toptable_clean.ALL_contrast.mark_seqc.gz
	ln -s $< $@

TD-edger.toptable.gz: ../../../RNAseq/dataset/v1/DGE_H9_TD/edger.toptable_clean.ALL_contrast.mark_seqc.gz
	link_install $< $@

GeneIDs_%-up: %-edger.toptable.gz
	bawk '$$6==1 {print $$2}' $< > $@
GeneIDs_NPC_TD-up: GeneIDs_NPC-up GeneIDs_TD-up
	cat $< $^2 | symbol_count | cut -f 1 > $@

GeneIDs_NPC_TD-up.tss2kb.bed: GeneIDs_NPC_TD-up $(TSS) $(CHROM_SIZE)
	filter_1col 4 $< < $^2 | bedtools slop -b 1000 -i stdin -g $^3 | bedtools sort -i stdin | bedtools merge -i stdin -c 4 -o distinct | bawk '{print $$1~3,$$4"_"i++}' > $@

%.bed.fa: $(BIOINFO_REFERENCE_ROOT)/gencode/dataset/hsapiens/32/GRCh38.primary_assembly.genome.clean_id.fa %.bed
	bedtools getfasta -name -fi $< -bed $^2 -fo $@


################################
#	genes association      #
################################

tss.slop: $(GENCODE_DIR)/primary_assembly.annotation.tss
	bedtools sort -i $< | bedtools slop -l $(REGREGION_DISTANCE) -r $(REGREGION_DISTANCE) -i stdin -g $(GENCODE_DIR)/chrom.info | bedtools sort -i stdin > $@

%.regulated_genes.bed: %.bed tss.slop
	bedtools intersect -a $< -b $^2 -loj > $@


################################
#     genhancer association    #
################################

%.tts_genom_coords.genehancer.gz: %.tts_genom_coords.bed $(BIOINFO_REFERENCE_ROOT)/GeneHancer/dataset/v5/genehancer.connected_genes.all.clean_ensg.rr_coords.tissue.bed
	bedtools intersect -a $< -b $^2 -wa -wb -loj | gzip > $@


hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.genehancer.bed: hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.bed /sto1/ref/bioinfotree/task/GeneHancer/dataset/v5/genehancer.connected_genes.all.clean_ensg.rr_coords.bed
	bedtools intersect -a $< -b $^2 -wa -wb -loj > $^3
stats.rr_type: hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.genehancer.gz
	bawk '{print $$1~3,$$7,$$21}' $< | uniq | cut -f 4,5 | symbol_count > $^2
stats.target_genes: hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.genehancer.gz
	bawk '{print $$7,$$23}' $< | sort | uniq | cut -f 1 | symbol_count > $^2
stats.n_rr: hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.genehancer.gz
	bawk '{print $$1~3,$$7}' $< | uniq | cut -f 4 | symbol_count > $^2



#############
#   grafo   #
#############

hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.best.tsv: hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.bed
	find_best 7:15 5 < $< | cut -f 5,7,15 > $@

.META: hg38-*.signature_only.bed.fa.tpx.tts_genom_coords.best.tsv
	1	score
	2	sequence_id
	3	region_name

%.header_added.gz: %.gz
	(bawk -M $< | cut -f 2 | transpose; zcat $< ) | gzip > $@

%.header_added: %
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@

hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.best.ccREtype.tsv: hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.best.tsv /sto1/ref/bioinfotree/task/encode-screen/dataset/v13/hg38-NPC.bed
	cat $< | translate -a -r <(cut -f4,10 $^2) 3 > $@

.META: hg38-*.signature_only.bed.fa.tpx.tts_genom_coords.best.ccREtype.tsv
	1	score
	2	sequence_id
	3	region_name
	4	cCRE_type

lncRNA_network.pre: hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.best.tsv
	cat $< | translate -a -d -j -f 3 $< 3 > $@

.META: lncRNA_network.pre
	1	score_lnc1_cCRE
	2	lnc1
	3	cCRE
	4	score_lnc2_cCRE
	5	lnc2

lncRNA_network.cutoff%.net: lncRNA_network.pre
	bawk -v C=$* '$$score_lnc1_cCRE>=C && $$score_lnc2_cCRE>=C {print $$lnc1,$$lnc2,$$cCRE}' $< | sort | uniq | cut -f 1,2 | symbol_count | bsort -k3,3nr > $@

lncRNA_network.point.cutoff%.net: lncRNA_network.pre
	bawk -v C=$* '($$score_lnc1_cCRE==C || $$score_lnc2_cCRE==C) && !($$score_lnc1_cCRE<C || $$score_lnc2_cCRE<C){print $$lnc1,$$lnc2,$$cCRE}' $< | sort | uniq | cut -f 1,2 | symbol_count | bsort -k3,3nr > $@

lncRNA_network.cutoff.matrix: lncRNA_network.cutoff15.net lncRNA_network.point.cutoff11.net lncRNA_network.point.cutoff12.net lncRNA_network.point.cutoff13.net lncRNA_network.point.cutoff14.net
	matrix_reduce -l '$^' 'lncRNA_network.*.net' | fasta2tab | bawk '$$2>=$$3 {print $$2" "$$3,$$1,$$4}' | sed 's/point.//; s/cutoff//' | tab2matrix -r lncRNAcouple > $@

######################
#
# da rivedere
#

net: hg38-NPC.signature_only.bed.fa.tpx.tts_genom_coords.best.tsv
	bawk '{print $$2,$$3,$$1}' $< > $@

net.jaccard: net
	./jaccard_bipartite.py < $< > $@

net.jaccard.plot: net.jaccard
	bawk '$$1>$$2 {print $$1" "$$2,$$3,$$4,$$5,$$6}' $< > $@
net.jaccard.grafo: net.jaccard
	bawk '$$3==15 && $$1>$$2 && $$6>0 {print $$1,$$2,$$6}' $< > $@

.META: net.jaccard.grafo
	1	seq1
	2	seq2
	3	jaccard

ssRNA.len: ssRNA.fa
	fastalen < $< > $@
