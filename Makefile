git_dir="/Users/mizuhirosuzuki/Documents/GitHub/revisit_crc"

simulation:
	Rscript R/CRC_no_covariate.R $(git_dir)
	Rscript R/CRC_cov.R $(git_dir)

write_up:
	cd Paper; \
		pdflatex paper.tex; \
		bibtex paper; \
		pdflatex paper.tex; \
		pdflatex paper.tex

all:
	make simulation
	make write_up
