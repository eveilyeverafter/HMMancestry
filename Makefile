all: install

test: install
	Rscript -e 'library(methods); devtools::test()'

document: roxygen staticdocs

roxygen:
	@mkdir -p man
	Rscript -e "library(methods); devtools::document()"

install:
	R CMD INSTALL --no-test-load .


build:
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr HMMancestry*gz | tail -n1`
	@rm -f `ls -1tr HMMancestry*gz | tail -n1`
	@rm -rf HMMancestry.Rcheck

doc: 
	if [ -e HMMancestry.pdf ]; then rm HMMancestry.pdf; fi
	R CMD Rd2pdf -M ./ --title='HMMancestry' -o "HMMancestry.pdf"
# No real targets!
.PHONY: all test document install