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
	R CMD check --no-manual `ls -1tr fbgenotyper*gz | tail -n1`
	@rm -f `ls -1tr fbgenotyper*gz | tail -n1`
	@rm -rf fbgenotyper.Rcheck

# No real targets!
.PHONY: all test document install