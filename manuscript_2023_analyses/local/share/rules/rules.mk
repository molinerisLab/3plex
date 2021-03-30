TRIPLEXATOR_PARAM=-l 10 -L -1 -e 20 -E -1

ssRNA.fa: /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/gencode.v32.transcripts.fa.gz ../../local/share/data/selected_lncRNA
	zcat $< | fasta2oneline | tr "|" "\t" | filter_1col 6 <(cut -f 1 $^2) | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,10 | tab2fasta | fold > $@
ssRNA.fa: /sto1/ref/bioinfotree/task/gencode/dataset/hsapiens/32/gencode.v32.transcripts.fa.gz ../../local/share/data/selected_lncRNA
	zcat $< | fasta2oneline | tr "|" "\t" | filter_1col 6 <(cut -f 1 $^2) | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,10 | tab2fasta | fold > $@

peak.fa: /sto1/epigen/Fatemeh_hESC/dataset/PARCLIP_P300_v1/Fatemeh-v2-P300_merged_reps.gencode_exon_annotation.seq.header_added.gz ../../local/share/data/selected_lncRNA
	bawk 'NR>1 {print $$exon_geneID, $$peak_seq}' $< | filter_1col 1 $^2 | id2count -b 1 | bawk '{print ">"$$1; print $$2}' > $@

%.fa.tfo: ssRNA.fa
	ssh epigen \
	sudo su - edoardo \
	/home/edoardo/src/triplexator/bin/triplexator $(TRIPLEXATOR_PARAM) -fm 0 -of 0 -o $@ -rm 2 -p 24 -ss $<

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
