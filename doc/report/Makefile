#****************************************************************************
#
# Makefile for pdf article.
#
#****************************************************************************

INPUT:=ds
EXTRA:=ds-non-paper

pdf:
	pdflatex -interaction=nonstopmode ${INPUT}.tex
	bibtex ${INPUT}.aux
	pdflatex -interaction=nonstopmode ${INPUT}.tex
	pdflatex -interaction=nonstopmode ${INPUT}.tex
	evince ${INPUT}.pdf

latex:
	latex -interaction=nonstopmode ${INPUT}.tex
	bibtex ${INPUT}
	latex -interaction=nonstopmode ${INPUT}.tex
	latex -interaction=nonstopmode ${INPUT}.tex
	dvipdfm ${INPUT}.dvi
	evince ${INPUT}.pdf

clean:
	- rm -f ${INPUT}.dvi *.aux *.bbl ${INPUT}.blg ${INPUT}.log ${INPUT}.out ${INPUT}.toc ${INPUT}.idx *.cpt *.aux *.lof *.lot *.brf *~

