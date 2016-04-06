from pylatex import Document, Section, Tabular, Package, Command, Figure
from pylatex.utils import NoEscape, bold
import os


class CreatePDF(object):

    def __init__(self, tilemetrics, controlmetrics, errormetrics, extractionmetrics, indexmetrics, qualitymetrics,
                 corintmetrics):
        self.tile = tilemetrics
        self.control = controlmetrics
        self.error = errormetrics
        self.extract = extractionmetrics
        self.index = indexmetrics
        self.quality = qualitymetrics
        self.corint = corintmetrics

    def create_pdf(self):
        pdflatex = '/usr/local/texlive/2015/bin/x86_64-linux/pdflatex'
        doc = Document()
        doc.packages.append(Package('geometry', options=['tmargin=0.75in', 'lmargin=0.75in', 'rmargin=0.75in']))
        doc.packages.append(Package('datetime', options=['ddmmyyyy']))
        doc.preamble.append(Command('newdateformat',
                                    NoEscape(r'mydate}{\twodigit{\THEDAY}/\twodigit{\THEMONTH}/\THEYEAR')))
        doc.preamble.append(Command('title', 'Quality Assessment'))
        doc.preamble.append(Command('date', NoEscape(r'\mydate\today')))
        doc.append(NoEscape(r'\maketitle'))
        doc.append(Command('begin', 'flushright'))
        doc.append(Command('Large', bold('<Run name>')))
        doc.append(Command('end', 'flushright'))

        avg_qual = self.get_avg_qual(self.quality)
        self.get_qual_graph(self.quality)

        with doc.create(Section('Quality data')):
            with doc.create(Tabular('l|l|l')) as table:
                table.add_row(('Data', 'Value', 'Pass/Fail'))
                table.add_hline()
                table.add_row(('Mean Cluster Density (k/mm2)', format(self.tile.mean_cluster_density/1000, '.2f'), ''))
                table.add_row(('Clusters passed filter (%)', '%s%%' % (format(self.tile.percent_pf_clusters, '.2f')),
                               ''))
                table.add_row(('Average >= Q30', '%s%%' % (format(avg_qual, '.2f')), ''))
                table.add_row(('1st full read >= Q30', '', ''))
                table.add_row(('2nd full read >= Q30', '', ''))
        Figure.add_plot()

        '''
        with doc.create(Section('Tile Metrics')):
            doc.append(self.tile)
        with doc.create(Section('Control Metrics')):
            doc.append(self.control)
        with doc.create(Section('Error Metrics')):
            doc.append(self.error)
        with doc.create(Section('Extraction Metrics')):
            doc.append(self.control)
        with doc.create(Section('Index Metrics')):
            doc.append(self.index)
        with doc.create(Section('Quality Metrics')):
            doc.append(self.quality)
        '''

        doc.generate_pdf('output', clean_tex=False, compiler=pdflatex)
        os.system("xdg-open /home/cuser/PycharmProjects/AMLpipeline/output.pdf")

    def get_avg_qual(self, qualitymetrics):
        col_list = list(qualitymetrics.df)
        col_list.remove('cycle')
        col_list.remove('lane')
        col_list.remove('tile')
        qualitymetrics.df['total'] = qualitymetrics.df.sum(axis=1)
        total_clusters = qualitymetrics.df['total'].sum(axis=0)

        col_list = range(29, 50)
        a = qualitymetrics.df[col_list].sum(axis=1)
        avg_qual = (float(sum(a)) / float(total_clusters))*100

        return avg_qual

    def get_qual_graph(self, qualitymetrics):
        quality_df = qualitymetrics.df

        cols_to_drop = ['cycle', 'lane', 'tile']
        quality_df = quality_df.drop(cols_to_drop, axis=1)
        qual_scores = quality_df.sum(axis=0).values
        qual_list = map(int, qual_scores.tolist())
        less_than_q30 = qual_list[0: 28]
        more_than_q30 = qual_list[29: 50]

        myint = 1000000
        y1 = [x / myint for x in less_than_q30]
        y2 = [x / myint for x in more_than_q30]
        x1 = range(1, 29)
        x2 = range(29, 50)

        plt.bar(x1, y1, color='b')
        plt.bar(x2, y2, color='g')
        plt.ylabel('Clusters (millions)')
        plt.xlabel('Q score')
        plt.show()