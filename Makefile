all: y-fuse.html

data:
	Rscript make/dat-fish.R make/dat-squa.R

fits: output/results/fish-6par.rds output/results/squa-6par.rds

output/results/fish-6par.rds: output/data/fish.rds
	Rscript R/fit-fish-6par.R

output/results/squa-6par.rds: output/data/squa.rds
	Rscript R/fit-squa-6par.R

y-fuse.R: output/results/fish-6par.rds output/results/squa-6par.rds

y-fuse.Rmd: y-fuse.R
	Rscript -e "library(sowsear); sowsear('$<', 'Rmd')"
y-fuse.md: y-fuse.Rmd
	Rscript -e "library(knitr); knit('$<')"
y-fuse.html: y-fuse.md
	Rscript -e "library(markdown);\
	 opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images');\
	 markdownToHTML('$<', '$@', options=opts)"

clean-data: 
	rm -f output/data/fish.rds
	rm -f output/data/squa.rds

clean-outputs:
	rm -rf output/figs
	rm -rf output/results

deps:
	Rscript make/dependencies.R
