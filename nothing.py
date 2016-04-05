from pylatex import Document, Section, Tabular, Package, Subsection, Math
from pylatex.utils import italic

doc = Document()
doc.packages.append(Package('geometry', options=['tmargin=1cm', 'lmargin=10cm']))

with doc.create(Section('The simple stuff')):
    doc.append('Some regular text and some')
    doc.append(italic('italic text. '))
    doc.append('\nAlso some crazy characters: $&#{}')
    with doc.create(Subsection('Math that is incorrect')):
        doc.append(Math(data=['2*3', '=', 9]))

    with doc.create(Subsection('Table of something')):
        with doc.create(Tabular('rc|cl')) as table:
            table.add_hline()
            table.add_row((1, 2, 3, 4))
            table.add_hline(1, 2)
            table.add_empty_row()
            table.add_row((4, 5, 6, 7))

doc.generate_pdf('full.pdf')
