merlo_genes:
	---------- Forwarded message ---------
	From: Carla Liaci <carla.liaci@unito.it>
	Date: Tue, Oct 12, 2021 at 11:07 AM
	Subject: lista per predizione triplex
	To: <ivan.molineris@unito.it>, Annalaura Tamburrini <annalaura.tamburrini@unito.it>
	Cc: Giorgio Roberto Merlo <giorgioroberto.merlo@unito.it>

merlo_genes.in_gencode: merlo_genes ../../../dataset/v1/longest_transcripts.fa
	cut -f 2 $< | unhead -n 2 | translate -a -v -e 0 <(grep '>' $^2 | tr -d '>' | append_each_row 1) 1 > $@

merlo_genes.clean: merlo_genes.in_gencode
	bawk '$$2==1 {print $$1}' $< > $@
