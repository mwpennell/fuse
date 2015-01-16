all: fusion.html

data:
	Rscript make/dat-fish.R make/dat-squa.R make/dat-fish-indchrom.R make/dat-squa-indchrom.R

fusion.R: output/results/chromosome-fish.rds output/results/chromosome-squa.rds output/results/karyotype-fish.rds output/results/karyotype-squa.rds output/results/karyotype-squa-6par.rds output/results/counts-raw.csv

fusion.Rmd: fusion.R
	Rscript -e "library(sowsear); sowsear('$<', 'Rmd')"
fusion.md: fusion.Rmd
	Rscript -e "library(knitr); knit('$<')"
fusion.html: fusion.md
	Rscript -e "library(markdown);\
	 opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images');\
	 markdownToHTML('$<', '$@', options=opts)"

clean-data: 
	rm -f output/data/fish.rds
	rm -f output/data/squa.rds
	rm -f output/data/fish-indchrom.rds
	rm -f output/data/squa-indchrom.rds

clean-outputs:
	rm -rf output/figs
	rm -rf output/results

deps:
	Rscript make/dependencies.R
