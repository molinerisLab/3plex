ssRNA2.fa: $(GENCODE_DIR)/gencode.v32.transcripts.fa.gz selected_lncRNA2
	zcat $< | fasta2oneline | tr "|" "\t" | filter_1col 6 <(cut -f 1 $^2) | bawk '$$8!="retained_intron"' | find_best 6 7 | cut -f 6,10 | tab2fasta | fold > $@

