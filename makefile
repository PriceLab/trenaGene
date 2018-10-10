all:  docs install

docs:
	R -e "devtools::document()"
vig:
	R -e "devtools::build_vignettes()"

build:
	(cd ..; R CMD build --no-build-vignettes trenaGene)

install:
	(cd ..; R CMD INSTALL trenaGene)

check:
	(cd ..; R CMD check `ls -t trenaGene_* | head -1`)

biocCheck:
	(cd ..; R CMD BiocCheck `ls -t trenaGene_* | head -1`)
