CONDA_ROOT?=/opt/conda
CONDA_VERSION?=miniconda3
CONDA_ACTIVATE=source $(CONDA_ROOT)/$(CONDA_VERSION)/etc/profile.d/conda.sh; conda activate

####################
#
# Generic rules
#

%.header_added: %
	(bawk -M $< | cut -f 2 | transpose; cat $< ) > $@
%.header_added.gz: %.gz
	(bawk -M $< | cut -f 2 | transpose; zcat $< ) | gzip > $@
%.xlsx: %.gz
	zcat $< | tab2xlsx > $@
%.xlsx: %
	cat $< | tab2xlsx > $@




####################
#
# Corrected rules
#

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.AUC_comp.ALL_versions_peralta.gz:
	@echo 'downloaded from peralta: /hpcnfs/data/epi/epigen/TPXcCRE/dataset/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.AUC_comp.ALL_versions_peralta.gz'

tpx_paramspace_AUC.method.Npeaks_version.technique.gz: v8.3_ReChIRP_overlap/tpx_paramspace_AUC.all_human_noSingleNt.method.Npeaks_version.gz v8.5_ReChIRP_overlap_mouse/tpx_paramspace_AUC.all_mouse_noSingleNt.method.Npeaks_version.gz /sto1/epigen/ReChIRP/ChIP_ENCODE_pipeline/dataset/ssRNA_technique version_method.dir
	zcat $< $^2 | translate -a $^3 2 | translate $^4 1 | gzip > $@

tpx_paramspace_AUC.method.Npeaks_version.technique.triplex_ssRNA.gz: tpx_paramspace_AUC.method.Npeaks_version.technique.gz triplex_ssRNA
	zcat $< | filter_1col 2 $^2 | gzip > $@

tpx_paramspace_AUC.method.Npeaks_version.technique.triplex_ssRNA.triplexator_default_params.gz: tpx_paramspace_AUC.method.Npeaks_version.technique.triplex_ssRNA.gz
	bawk '$$5 == "ss0" && $$7 == 20 && $$8 == 40 && $$9 == "on" && $$10 == 1 && $$11 == "t_pot_norm"' $< | gzip > $@

.META: tpx_paramspace_AUC.method.Npeaks_version*.gz
	1	peak_method_all	idr_conservative
	2	ssRNA	7SK
	3	technique	ChIRP-seq
	4	n_peaks	249
	5	singleStrandedness	ss0
	6	min_len	10
	7	error_rate	10
	8	guanine_rate	10
	9	repeat_filter	off
	10	consecutive_error	1
	11	predictor	Score_best
	12	AUC	0.537313914291705

best_single_params.tsv: tpx_paramspace_AUC.method.Npeaks_version.technique.gz
	zcat $< | find_best 2 12 > $@


####################
#
# Anova Exclude
#

%.anova_exclude.gz: %.gz
	export OPENBLAS_NUM_THREADS=1; seq 1 10000 | parallel -j 50 --group ./anova_perm.sh > $@.tmp; \
        grep -v '^>' $@.tmp | gzip > $@; \
        rm $@.tmp

#FORMULA=AUC~peak_method_all+ssRNA+n_peaks+min_len
FORMULA=AUC~peak_method_all+technique+n_peaks+singleStrandedness+min_len+error_rate+guanine_rate+repeat_filter+consecutive_error+predictor

%.true_prediction: %.gz
	bmodel -i -a $< $(FORMULA) | grep -v '^>' > $@




####################
#
# ROC curves
#

best_single_params.triplex_ssRNA_scores.peralta.gz:
	@echo 'downloaded from peralta /hpcnfs/data/epi/epigen/TPXcCRE/dataset'

best_single_params.triplex_ssRNA_scores.gz: best_single_params.triplex_ssRNA_scores.peralta.gz
	zcat $< | grep -v Duplex_ID | bawk '{split($$1,a,";"); print a[2],$$0}' | cut -f1,3- | gzip > $@

.META: best_single_params.triplex_ssRNA_scores.gz
	1	ssRNA
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



tpx_paramspace_AUC_cmp.human_mouse.gz: v8.6_ReChIRP_idr_overlap_top1000/tpx_paramspace_AUC_cmp.gz v8.7_ReChIRP_idr_overlap_top1000_mouse/tpx_paramspace_AUC_cmp.gz
	cat $< $^2 > $^@

tpx_paramspace_AUC_cmp.triplex_ssRNA.mean_AUC.gz: tpx_paramspace_AUC_cmp.human_mouse.gz selected_ssRNA.triplex_ssRNA
	zgrep -v TTS_covered_frac $< | bawk 'BEGIN{FS = "\t";OFS = ";"}{print $$1"\t"$$2,$$4~10"\t"$$12}' | sort -u | filter_1col 1 $^2 | cut -f2,3 | sort -k1,1 | stat_base -g -a | sort -k2,2nr | gzip > $^3

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.human_mouse.gz: v8.6_ReChIRP_idr_overlap_top1000/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.gz v8.7_ReChIRP_idr_overlap_top1000_mouse/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.gz
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
	zcat $< | ../local/src/ROC.R neg_pos $$(zcat $< | head -n1 | cut -f2- | tr "\t" " ") -d "<" | gzip > $@

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.best_param_setting.gz:
	matrix_reduce -t 'raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.*.AUC_cmp.gz' | grep -v 'pred_1' | bawk '{print $$1,$$2,$$4}' | sort -k3,3nr | uniq | bawk '{split($$2,a,"_"); print $$1,"ss"a[2],a[3],a[5],a[6],a[7],a[8],$$3}' | gzip > $@

raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.Stability_norm_undercount.best_param_setting.gz: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.ALL_ssRNA.human_mouse.gz selected_ssRNA.triplex_ssRNA
	bawk '{split($$1,a,";"); print a[1],a[2]"_"a[4]"_1_"a[6]"_"a[7]"_"a[8]"_"a[9], $$2,$$4, $$Stability_norm_undercount}' $< | filter_1col 1 $^2 | grep -w 0_8_1_20_70_off_3 | bawk 'BEGIN{print "peak","neg_pos","0_8_1_20_70_off_3"}{print $$3,$$4,$$5}' | gzip > $@

v8.6_v8.7_method_cmp.matrix.gz: v8.6_ReChIRP_idr_overlap_top1000/ROC/tpx.neg_pos.ALL_scores.gz v8.7_ReChIRP_idr_overlap_top1000_mouse/ROC/tpx.neg_pos.ALL_scores.gz raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.triplex_ssRNA.Stability_norm_undercount.best_param_setting.gz
	cat <(zcat $< | cut -f2,5,6 | sed 's/peak;ssRNA/peak/') <(zcat $^2 | bawk 'NR>1{print $$2,$$5,$$6}') | translate -a <(zcat $^3) 1 | gzip > $@

v8.6_v8.7_method_cmp.matrix.AUC_cmp.tsv: v8.6_v8.7_method_cmp.matrix.gz
	zcat $< | ../local/src/ROC.R neg_pos $$(zcat $< | head -n1 | cut -f3- | tr "\t" " ") -d "<" > $@



####################
#
# TO FIX
#

tpx_paramspace_AUC.idr_overlap_top1000.human.best_general_params.gz: tpx_paramspace_AUC.idr_overlap_top1000.human.gz best_general_params.tsv
	bawk '{print $$2";"$$3";"$$4";"$$5";"$$6";"$$7, $$0}' $< | filter_1col 1 $^2 | cut -f2- | gzip > $^3

best_single_params.ALL_ssRNA_scores.uv2000.gz: v8_ChIRP_neg_rand/tpx_paramspace/SRA1_ss0_unpairedWindow/SRA1.neg_pos_rand.bed/min_length~8/max_length~-1/error_rate~20/guanine_rate~40/filter_repeat~off/consecutive_errors~1/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic
	matrix_reduce -t -l 'v8_ChIRP_neg_rand/tpx_paramspace/MEG3_ss50_unpairedWindow/MEG3.neg_pos_rand.bed/min_length~10/max_length~-1/error_rate~20/guanine_rate~40/filter_repeat~off/consecutive_errors~3/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic $<' '*/tpx_paramspace/*_*_unpairedWindow/*.neg_pos_rand.bed/min_length~*/max_length~-1/error_rate~*/guanine_rate~*/filter_repeat~*/consecutive_errors~*/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic' | gzip > $@
best_single_params.ALL_ssRNA_scores.gz: best_single_params.ALL_ssRNA_scores.peralta.gz best_single_params.ALL_ssRNA_scores.uv2000.gz
	zcat $< $^2 | grep -v Duplex_ID | bawk '{split($$1,a,";"); print a[2],$$0}' | cut -f1,3- | gzip > $@


.META: best_single_params.ALL_ssRNA_scores.gz
	1	ssRNA
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
	16	Stability_norm_overcount
	17	Stability_norm_undercount
	18	PROB__fitted_model
best_single_params.ALL_ssRNA_scores.human.header_added.gz: best_single_params.ALL_ssRNA_scores.header_added.gz v8.3_ReChIRP_overlap/human_lncRNA
	zcat $< | filter_1col --header 1 1 $^2 | gzip > $^3
best_single_params.ALL_ssRNA_scores.mouse.header_added.gz: best_single_params.ALL_ssRNA_scores.header_added.gz v8.5_ReChIRP_overlap_mouse/mouse_lncRNA
	zcat $< | filter_1col --header 1 1 $^2 | gzip > $^3
best_single_params.ALL_ssRNA_scores.triplex_ssRNA.header_added.gz: best_single_params.ALL_ssRNA_scores.header_added.gz triplex_ssRNA
	zcat $< | filter_1col --header 1 1 $^2 | gzip > $^3
tpx_paramspace_AUC.idr_conservative.best_general_params.triplex_ssRNA.header_added.gz: tpx_paramspace_AUC.idr_conservative.mouse.best_general_params.header_added.gz tpx_paramspace_AUC.idr_conservative.human.best_general_params.header_added.gz triplex_ssRNA
	cat <(zcat $<) <(zcat $^2 | unhead) | filter_1col --header 1 1 $^3 | gzip > $^4

BIN=20
avg_roc_pre.bin: avg_roc_pre.tsv
	bawk 'NR>1 {print $$2}' $< | round_table -p 4 | sort -k1,1g | uniq | enumerate_rows | perl -lane 'BEGIN{$$,="\t"; $$i=0; print ">$$i"; $$i++} $$F[0] = $$F[0] % int((2216/$(BIN))-1); if($$F[0]==0){print ">$$i"; $$i++} print $$F[1]' | fasta2tab > $@
avg_roc_pre.bin.avg: avg_roc_pre.bin
	stat_base -o -g -a < $< > $@
avg_roc_pre.tsv.avg1: avg_roc_pre.tsv avg_roc_pre.bin avg_roc_pre.bin.avg
	round_table -p 4 $< | uniq | translate -a -f 2 <(echo -e "spec_bin\tspec"; cat $^2) 2 | translate -a <(echo -e "spec_bin\tspec_bin_avg"; cat $^3) 3 > $@

tpx_paramspace_AUC.method.Npeaks_version.technique.triplex_ssRNA.len8_10.gz: tpx_paramspace_AUC.method.Npeaks_version.technique.triplex_ssRNA.gz
	bawk '$$6 != 16 && $$6 != 12' $< | gzip > $@

best_single_params.triplex_ssRNA.tsv: tpx_paramspace_AUC.method.Npeaks_version.technique.gz triplex_ssRNA
	zcat $< | filter_1col 2 $^2 | grep -v TTS_covered_frac | find_best 2 12 > $@

best_single_params.triplex_ssRNA.AUC_comp.tsv: raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.AUC_comp.ALL_versions_peralta.gz filter_best_single_params version_method.dir
	bawk '{print $$1";"$$2";"$$3";"$$4";"$$5";"$$6";"$$7";"$$8";"$$9,$$0; print $$1";"$$2";"$$3";"$$4";"$$5";"$$6";"$$7";"$$8";"$$10,$$0}' $< | filter_1col 1 $^2 | grep -v 'TTS_covered_frac' | translate $^3 2 > $@
filter_best_single_params: best_single_params.triplex_ssRNA.tsv version_method.dir
	cat $< | translate -a -d -j -f 2 $^2 1 | cut -f2,3,6- | cut -f-9 | tr "\t" ";" > $@

filter_best_single_params.file_name.tsv: filter_best_single_params
	cat $< | tr ";" "\t" | cut -f-8 | bawk '{print $$1"/tpx_paramspace/"$$2"_"$$3"_unpairedWindow/"$$2".neg_pos_rand.bed/min_length~"$$4"/max_length~-1/error_rate~"$$5"/guanine_rate~"$$6"/filter_repeat~"$$7"/consecutive_errors~"$$8"/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.gz"}' > $@

best_single_params.no_triplex_ssRNA.tsv: triplex_ssRNA best_single_params.tsv
	filter_1col -v 2 $< < $^2 | bsort -k12,12nr > $@
best_single_params.ALL_ssRNA.tsv: tpx_paramspace_AUC.method.Npeaks_version.technique.header_added.gz triplex_ssRNA
	zcat $< | zgrep -v TTS_covered_frac | find_best -H 2 12 | sort -k12,12n | translate -a -v -e 0.5 <(bawk '{print $$1, 1}' $^2) 2 > $^3

tpx_paramspace_AUC.method.Npeaks_version.technique.best_ss.matrix: tpx_paramspace_AUC.method.Npeaks_version.technique.gz
	bawk '{print $$ssRNA,$$singleStrandedness,$$AUC}' $< | find_best 1:2 3 | tab2matrix -r lncRNA > $@

matrix.auc_no_cmp.method_species.all_versions.gz:
	zcat */*auc_no_cmp.method_species.gz | bawk 'NR==1 || $$1!="shuffling"' | gzip > $@
matrix.auc_no_cmp.method_species.all_versions.clean.gz: matrix.auc_no_cmp.method_species.all_versions.gz
	zcat $< | sed -e 's/triplexAligner/triplexAligner_NA_NA_NA_NA_NA_NA_NA_NA/g' -e 's/fasimLongtarget/fasim_LongTarget_NA_NA_NA_NA_NA_NA_NA_NA/g' | bawk '{split($$2,a,"_"); print $$1,a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[9],$$3,$$4}' | gzip > $^2








############################
# Genomic Peaks subsampling

v_test_randomiz_sub.%.summary.clean.gz:
	matrix_reduce -t 'v_test_randomiz_sub_*/*.$*.summary.clean.gz' | bawk 'NR==1{print $$2,$$3,$$4}NR>1 && $$2!="pos_neg" {print$$2,$$3,$$4}' | gzip > $@
v_test_randomiz_sub.summary.clean.tpx_method.gz:
	matrix_reduce -t 'v_test_randomiz_sub.*.summary.clean.gz' | grep -v pos_neg | cut -f 1,3,5 | bawk 'BEGIN{print "tpx_method","pos_neg","Stability_best"}{print}' | gzip > $@

AUC_singles.tsv: $(addsuffix .3plex.summary.clean.AUC_cmp.tsv, $(SAMPLES)) $(addsuffix .triplexAligner.summary.clean.AUC_cmp.tsv, $(SAMPLES)) $(addsuffix .fasimLongtarget.summary.clean.AUC_cmp.tsv, $(SAMPLES))
	matrix_reduce -t -l '$^' '*\.*\.summary.clean.AUC_cmp.tsv' | tr ";" "\t" | cut -f1,2,4,6 | grep -v AUC > $@
AUC_singles.matrix.xlsx: AUC_singles.tsv
	cat $< | cut -f1,2,4 | tab2matrix -r ssRNA | tab2xlsx > $@


##########
# Avg ROC
BIN=20
NEGATIVE_AMPLIFICATION=50
#NEGATIVE_AMPLIFICATION=1000
%.summary.clean.sens_spec.tsv: %.summary.clean.gz
	$(CONDA_ACTIVATE) /home/cciccone/.conda/envs/pROC_Env; \
	zcat $< | ../local/src/avgROC.R pos_neg pred2 --negative_subsampling $(NEGATIVE_AMPLIFICATION) > $@
v_test_randomiz_sub.summary.clean.sens_spec.tpx_method.gz:
	matrix_reduce -t 'v_test_randomiz_sub.*.summary.clean.sens_spec.tsv' | grep -v specificities | bawk 'BEGIN{print "tpx_method","iteration","specificities","sensitivities"}{print}' | gzip > $@
