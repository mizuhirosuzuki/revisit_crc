simulation:
	Rscript R/CRC_no_covariate.R
	Rscript R/CRC_cov.R
		$(data_dir) $(git_dir)

write_up:
	cd Paper; \
		pdflatex paper.tex; \
		bibtex paper; \
		pdflatex paper.tex; \
		pdflatex paper.tex

all:
	make simulation
	make write_up
