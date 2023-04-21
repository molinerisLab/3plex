CONDA_ROOT?=/opt/conda
CONDA_VERSION?=miniconda3
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate


#################
#  Generic rules

%.header_added: %
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@
%.header_added.gz: %.gz
	(bawk -M $< | cut -f 2 | transpose; zcat $< ) | gzip > $@
%.xlsx: %.gz
	zcat $< | tab2xlsx > $@
%.xlsx: %
	cat $< | tab2xlsx > $@



###############
#  ROC curves

triplexAligner.matrix.human_mouse.gz: $(addsuffix /triplexAligner.matrix.gz, $(V8_DIRS))
	r --vanilla -e 'h <- read.delim("$<",header=T);\\
m <- read.delim("$^2",header=T);\\
hm <- rbind(h,m[,colnames(h)]);\\
write.table(hm,file=gzfile("$@"),quote=F,sep="\t",row.name=F,col.names=T)'

fasim_LongTarget.matrix.human_mouse.gz: $(addsuffix /fasim_LongTarget.matrix.gz, $(V8_DIRS))
	r --vanilla -e 'h <- read.delim("$<",header=T);\\
m <- read.delim("$^2",header=T);\\
hm <- rbind(h,m[,colnames(h)]);\\
write.table(hm,file=gzfile("$@"),quote=F,sep="\t",row.name=F,col.names=T)'

d3plex.matrix.human_mouse.gz: $(addsuffix /d3plex.matrix.gz, $(V8_DIRS))
	r --vanilla -e 'h <- read.delim("$<",header=T);\\
m <- read.delim("$^2",header=T);\\
hm <- rbind(h,m[,colnames(h)]);\\
write.table(hm,file=gzfile("$@"),quote=F,sep="\t",row.name=F,col.names=T)'

%.matrix.human_mouse.auc_no_cmp.gz: %.matrix.human_mouse.gz
	set +u; source /opt/conda/miniconda3/etc/profile.d/conda.sh; conda activate ; conda activate  pROC_Env;\
	zcat $< | ../../local/src/ROC_no_comp.R pos_neg $$(zcat $< | head -n1 | cut -f2- | tr "\t" " ") -d "<" | gzip > $@

%.matrix.human_mouse.auc_no_cmp.method_species.gz: %.matrix.human_mouse.auc_no_cmp.gz
	bawk 'BEGIN{print"shuffling","tpx_method","auc","species"}{print $$1,"$*",$$2,"human_mouse"}' $< | gzip > $@

ifeq ($(SHUFFLE_METHOD),transcript)
d3plex.matrix.human_mouse.auc_no_cmp.method_species.gz: d3plex.matrix.human_mouse.auc_no_cmp.gz
	bawk 'BEGIN{print "shuffling","tpx_method","auc","species"}{split($$1,a,"_"); print a[1],"d3plex_"a[2]"_"a[3]"_"a[4]"_"a[5]"_"a[6]"_"a[7]"_"a[8]"_"a[9],$$2,"human_mouse"}' $< | gzip > $@
else
d3plex.matrix.human_mouse.auc_no_cmp.method_species.gz: d3plex.matrix.human_mouse.auc_no_cmp.gz
		bawk 'BEGIN{print"shuffling","tpx_method","auc","species"}{print "genomic_regions",$$1,$$2,"human_mouse"}' $< | gzip > $@
endif







tpx_paramspace_AUC_cmp.human_mouse.gz: $(addsuffix /tpx_paramspace_AUC_cmp.gz, $(V8_DIRS))
	cat $< $^2 > $@

tpx_paramspace_AUC_cmp.triplex_ssRNA.mean_AUC.gz: tpx_paramspace_AUC_cmp.human_mouse.gz selected_ssRNA.triplex_ssRNA
	zgrep -v TTS_covered_frac $< | bawk 'BEGIN{FS = "\t";OFS = ";"}{print $$1"\t"$$2,$$4~10"\t"$$12}' | sort -u | filter_1col 1 $^2 | cut -f2,3 | sort -k1,1 | stat_base -g -a | sort -k2,2nr | gzip > $^3

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.human_mouse.gz: $(addsuffix /raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.gz, $(V8_DIRS))
	cat $< $^2 > $@

META: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.human_mouse.gz 
	1	ssRNA;singleStrandedness;ssRNA;min_len;max_len;error_rate;guanine_rate;repeat_filter;consecutive_error
	2	Duplex_ID
	3	tpx_count_custom
	4	neg_pos
	5	tpx_count_standard
	6	t_pot_norm
	7	Duplex_length
	8	oligolength
	9	triplexator_norm_factor
	10	t_pot_custom
	11	TTS_covered_len
	12	TTS_covered_frac
	13	Stability_best
	14	Stability_tot_overcount
	15	Stability_tot_undercount
	16	Score_best
	17	Stability_norm_overcount
	18	Stability_norm_undercount


raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.%.matrix.gz: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.human_mouse.header_added.gz selected_ssRNA.triplex_ssRNA
	bawk '{split($$1,a,";"); print a[1],a[2]";"a[4]";"a[5]";"a[6]";"a[7]";"a[8]";"a[9], $$2";"$$4, $$$*}' $< | \
	filter_1col 1 $^2 | \
	bawk '{print $$3,"r_"$$2,$$4}' | \
	tab2matrix -r 'peak;ssRNA;neg_pos' | \
	bawk '{split($$1,a,";"); print a[3],$$0}' | \
	cut -f 1,3- | sed -e '1 s/;/_/g' -e '1 s/-//g' | \
	gzip > $@

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.%.AUC_cmp.gz: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.%.matrix.gz
	$(CONDA_ACTIVATE) pROC_Env; \
	zcat $< | ../../local/src/ROC.R neg_pos $$(zcat $< | head -n1 | cut -f2- | tr "\t" " ") -d "<" | gzip > $@

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.ALL_AUC.gz:
	matrix_reduce -t 'raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.*.AUC_cmp.gz' | grep -v 'pred_1' | bawk '{print $$1,$$2,$$4}' | sort -k3,3nr | uniq | gzip > $@

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.Stability_norm_undercount.best_param_setting.gz: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.human_mouse.gz selected_ssRNA.triplex_ssRNA
	bawk '{split($$1,a,";"); print a[1],a[2]"_"a[4]"_1_"a[6]"_"a[7]"_"a[8]"_"a[9], $$2,$$4, $$Stability_norm_undercount}' $< | filter_1col 1 $^2 | grep -w 0_8_1_20_70_off_3 | bawk 'BEGIN{print "peak","neg_pos","0_8_1_20_70_off_3"}{print $$3,$$4,$$5}' | gzip > $@

method_cmp.matrix.gz: $(addsuffix /ROC/tpx.neg_pos.ALL_scores.gz, $(V8_DIRS)) raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.human_mouse.gz selected_ssRNA.triplex_ssRNA
	cat <(zcat $< | cut -f2,4,5 | sed 's/peak;ssRNA/peak/') <(zcat $^2 | bawk 'NR>1{print $$2,$$4,$$5}') | translate -a <(bawk '{split($$1,a,";"); print a[1],a[2]"_"a[4]"_1_"a[6]"_"a[7]"_"a[8]"_"a[9], $$2,$$4, $$Stability_best}' $^3 | filter_1col 1 $^4 | grep -w 10_10_1_20_10_off_3 | bawk 'BEGIN{print "peak","neg_pos","d_3plex"}{print $$3,$$4,$$5}') 1 | gzip > $@

method_cmp.matrix.AUC_cmp.tsv: method_cmp.matrix.gz
	$(CONDA_ACTIVATE) pROC_Env;\
	zcat $< | ../../local/src/ROC.R neg_pos $$(zcat $< | head -n1 | cut -f3- | tr "\t" " ") -d "<" > $@

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.human_mouse.header_added.matrix.gz: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.human_mouse.header_added.gz selected_ssRNA.triplex_ssRNA
	bawk 'NR>1{split($$1,a,";"); id=a[4]"_"a[6]"_"a[7]"_"a[8]"_"a[9]"_"a[2]; print a[1],$$4";"$$2,id"_Stability_best",$$13; print a[1],$$4";"$$2,id"_Stability_norm",$$18}' $< | filter_1col 1 $^2 | cut -f2- | tab2matrix -r "pos_neg;peak;ssRNA" | bawk '{split($$1,a,";"); print a[1],$$0}' | cut -f 1,3- | gzip > $@

%.human_mouse.matrix.gz: ../v8.16.1_top1000_human_transcript_shuffle/%.matrix.gz ../v8.17.1_top1000_mouse_transcript_shuffle_k2/%.matrix.gz 
	r --vanilla -e 'h <- read.delim("$<",header=T);\\
m <- read.delim("$^2",header=T);\\
hm <- rbind(h,m[,colnames(h)]);\\
write.table(hm,file=gzfile("$@"),quote=F,sep="\t",row.name=F,col.names=T)'

triplexAligner.best_logE.pos_neg.ALL.human_mouse.matrix.auc_no_cmp.method_species.gz: triplexAligner.best_logE.pos_neg.ALL.human_mouse.matrix.auc_no_cmp.gz
	bawk 'BEGIN{print "shuffling","tpx_method","auc","species"}{print $$1,"triplexAligner",$$2,"human_mouse"}' $< | gzip > $@
fasim_LongTarget.best_meanStab.pos_neg.ALL.human_mouse.matrix.auc_no_cmp.method_species.gz: fasim_LongTarget.best_meanStab.pos_neg.ALL.human_mouse.matrix.auc_no_cmp.gz
	bawk 'BEGIN{print"shuffling","tpx_method","auc","species"}{print $$1,"fasimLongTarget",$$2,"human_mouse"}' $< | gzip > $@

