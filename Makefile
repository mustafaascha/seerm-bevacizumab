all: claims munge analysis bevpts results

claims: img dx mb str nsrg rad

munge: \
	cancers_loaded.rds\
	cancers_dx-img.csv.gz\
	cancers_nsrg.csv.gz\
	cancers_rad.csv.gz\
	cancers_joined-vars.csv.gz\
	cancers_prerecode.csv.gz\
	cancers_postrecode.csv.gz\
	cancers_w_mabs.csv.gz\
	cancers_before_exclusion.csv.gz\
	cancers.csv.gz\
  chemotherapy.csv.gz

R_OPTS=--vanilla

#am I using this correctly?
VPATH = \
	cache\
	cache/dx-imaging\
	cache/diagnoses\
	cache/neurosurg\
	cache/wbrt\
	cache/stereo\
	cache/mabs\
  cache/bevpts\
	extraction\
	munge\
	reports

# TODO: Figure out how to use eval, macros
# Functions
extract  = $(wildcard extraction/*$(1)*.R)
target   = $(subst extraction/,cache/$(2),$(subst R,csv.gz,$(1)))
targets  = $(foreach a,$(1),$(call target,$(a),$(2)/))
resource = extraction/$(subst csv.gz,R,$(notdir $1))

# implementation
img_srcs = $(call extract,img)
img_targets = $(call targets,$(img_srcs),dx-imaging)
img : $(img_targets)
$(img_targets) : $(img_srcs)
	Rscript $(call resource,$@)

dx_srcs = $(call extract,dx)
dx_targets = $(call targets,$(dx_srcs),diagnoses)
dx : $(dx_targets)
$(dx_targets) : $(dx_srcs)
	Rscript $(call resource,$@)

mb_srcs = $(call extract,mb)
mb_targets = $(call targets,$(mb_srcs),mabs)
mb : $(mb_targets)
$(mb_targets) : $(mb_srcs)
	Rscript $(call resource,$@)

str_srcs = $(call extract,str)
str_targets = $(call targets,$(str_srcs),stereo)
str : $(str_targets)
$(str_targets) : $(str_srcs)
	Rscript $(call resource,$@)

rad_srcs = $(call extract,rad)
rad_targets = $(call targets,$(rad_srcs),rad)
rad : $(rad_targets)
$(rad_targets) : $(rad_srcs)
	Rscript $(call resource,$@)

nsrg_srcs = $(call extract,nsrg)
nsrg_targets = $(call targets,$(nsrg_srcs),neurosurg)
nsrg : $(nsrg_targets)
$(nsrg_targets) : $(nsrg_srcs)
	Rscript $(call resource,$@)

# munging and postprocessing
cancers_loaded.rds:              load-data.R
	Rscript $<
cancers_dx-img.csv.gz:           claims-vars.R cancers_loaded.rds
	Rscript $<
cancers_nsrg.csv.gz:             treatment-nsrg.R cancers_dx-img.csv.gz
	Rscript $<
cancers_rad.csv.gz:              treatment-rad.R cancers_nsrg.csv.gz
	Rscript $<
cancers_joined-vars.csv.gz:      treatment-ster.R cancers_rad.csv.gz
	Rscript $<
cancers_prerecode.csv.gz:        claims-dates.R cancers_joined-vars.csv.gz
	Rscript $<
cancers_postrecode.csv.gz:       recoding.R cancers_prerecode.csv.gz
	Rscript $<
cancers_w_mabs.csv.gz:           mabs.R cancers_postrecode.csv.gz
	Rscript $<
cancers_before_exclusion.csv.gz: exclusion-vars.R cancers_w_mabs.csv.gz
	Rscript $<
cancers.csv.gz:                  misc.R cancers_before_exclusion.csv.gz
	Rscript $<
#exclusion.rds:                   exclusion.R cancers.csv.gz
#	Rscript $<

# retrieving information about lung cancer patients' treatment
bvsource = $(subst bk.,,extraction/$(subst csv.gz,R,$(notdir $1)))
bv_srcs = $(wildcard extraction/bev-patients*.R)
bv_targets = $(subst cache//,cache/bevpts/,$(call targets,$(bv_srcs),))
$(bv_targets) : $(bv_srcs) cancers.csv.gz
	Rscript $(call bvsource,$@)
split: split-bev.sh $(bv_targets)
	bash $< 
bevpts : $(bv_targets) split

# postprocessing, producing analyzable data
cancers-analytic.csv.gz : load_exclude.R bevpts
	Rscript $<
chemotherapy.csv.gz: chemo-courses.R bevpts
	Rscript $<
chemo-analytic.rds: chemo-analysis.R chemotherapy.csv.gz
	Rscript $< 

analysis: analysis.R bevpts cancers-analytic.csv.gz chemo-analytic.rds
	Rscript $< 
 
bev-results.pdf : bev-results.Rmd analysis chemo-analytic.rds exclusion.rds
	Rscript -e 'rmarkdown::render("bev-results.Rmd", output_format = "word_document")'

results: bev-results.pdf

clean: rm -f *.aux *.bbl *.blg *.log *.bak *~ *.Rout */*~ */*.Rout */*.aux */*.log

debug: 
	@echo bv targets: $(bv_targets)
	@echo
	@echo bv srcs: $(bv_srcs)
	@echo
