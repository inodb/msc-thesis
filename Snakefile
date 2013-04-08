DOCUMENTS = ['report']
TEXS = [doc+".tex" for doc in DOCUMENTS]
PDFS = [doc+".pdf" for doc in DOCUMENTS]
#FIGURES = ['fig1.pdf']

include: 'tex.smrules'

rule all:
    input: PDFS

rule zipit:
    output: 'upload.zip'
    input: TEXS, PDFS #, FIGURES
    shell: 'zip -T {output} {input}'

rule pdfclean:
    shell: "rm -f  {PDFS}"
