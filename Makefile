ps:
	fixrefurl.pl ../ref.bib
	fixbib.pl ../ref.bib
	latex ms
	bibtex ms
	latex ms
	latex ms
	dvips ms.dvi
	ps2pdf ms.ps
	open ms.pdf

qps:
	latex ms
	bibtex ms
	latex ms
	dvips ms.dvi
	ps2pdf ms.ps
	open ms.pdf

nobib:
	latex ms
	latex ms
	latex ms
	dvips ms.dvi
	ps2pdf ms.ps
	open ms.pdf

