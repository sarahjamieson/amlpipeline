from pylatex import Document, Section, Tabular, Package, Command, Figure, SubFigure, Subsection
from pylatex.utils import NoEscape, bold
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
        doc.packages.append(Package('needspace'))
        doc.preamble.append(Command('newdateformat',
                                    NoEscape(r'mydate}{\twodigit{\THEDAY}/\twodigit{\THEMONTH}/\THEYEAR')))
        doc.preamble.append(Command('title', 'MiSeq quality checks'))
        doc.preamble.append(Command('date', NoEscape(r'\mydate\today')))
        doc.append(NoEscape(r'\maketitle'))
        doc.append(Command('begin', 'flushright'))
        doc.append(Command('Large', bold('<Run name>')))
        doc.append(Command('end', 'flushright'))

        avg_qual = self.get_avg_qual(self.quality)
        self.get_qual_graph(self.quality, avg_qual)

        with doc.create(Section('Quality data')):
            with doc.create(Tabular(NoEscape(r'p{5cm}|c|c'))) as table:
                table.add_row(('Data', 'Value', 'Pass/Fail'))
                table.add_hline()
                table.add_row(
                    ('Mean Cluster Density (k/mm2)', format(self.tile.mean_cluster_density / 1000, '.2f'), ''))
                table.add_row(('Clusters passed filter (%)', '%s%%' % (format(self.tile.percent_pf_clusters, '.2f')),
                               ''))
                table.add_row(('Average >= Q30', '%s%%' % (format(avg_qual, '.2f')), ''))
                table.add_row(('1st full read >= Q30', '', ''))
                table.add_row(('2nd full read >= Q30', '', ''))

            with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))):
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_plot()
                    plot.add_caption('Q-score distribution plot (all reads all cycles)')
                self.get_qual_heatmap(self.quality)
                with doc.create(SubFigure()) as plot:
                    plot.add_plot()
                    plot.add_caption('Q-score heat map')
        with doc.create(Section('Phas/Prephas data')):
            with doc.create(Tabular(NoEscape(r'p{5cm}|c|c'))) as table:
                table.add_row(('Data', 'Value', 'Pass/Fail'))
                table.add_hline()
                table.add_row(('1st full read', '', ''))
                table.add_row(('2nd full read', '', ''))
        sample_id, index1, index2, percent_clusters, percent_pf, total_aligned_clusters, pf_aligned_clusters = \
            self.get_clusters_per_sample(self.index, self.tile)
        total_samples = len(sample_id)
        doc.append(Command('needspace', '10em'))
        with doc.create(Section('Indexing')):
            with doc.create(Subsection('Reads mapped to Index')):
                with doc.create(Tabular(NoEscape(r'c|c|c|c|c'))) as table:
                    table.add_row(('Total Reads', 'PF Reads', '% Reads Identified (PF)', 'Min', 'Max'))
                    table.add_hline()
                    table.add_row(
                        ('%s' % int(total_aligned_clusters), '%s' % int(pf_aligned_clusters), '%s' % percent_pf,
                         '%s' % min(percent_clusters), '%s' % max(percent_clusters)))
            with doc.create(Subsection('Reads per sample')):
                with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))):
                    with doc.create(Tabular(NoEscape(r'c|c|c|c'))) as table:
                        table.add_row(('Sample_ID', 'Index', 'Index2', '% Reads  Identified (PF)'))
                        table.add_hline()
                        item = 0
                        while item < total_samples:
                            table.add_row(('%s' % sample_id[item], '%s' % index1[item], '%s' % index2[item], '%s'
                                           % percent_clusters[item]))
                            item += 1
                    with doc.create(SubFigure()) as plot:
                        plot.add_plot()

        doc.generate_pdf('output', clean_tex=False, compiler=pdflatex)
        # os.system("xdg-open /home/cuser/PycharmProjects/AMLpipeline/output.pdf")

    def get_avg_qual(self, qualitymetrics):
        quality_df = qualitymetrics.df
        col_list = list(quality_df)
        col_list.remove('cycle')
        col_list.remove('lane')
        col_list.remove('tile')
        quality_df['total'] = quality_df.sum(axis=1)
        total_clusters = quality_df['total'].sum(axis=0)

        col_list = range(29, 50)
        a = quality_df[col_list].sum(axis=1)
        avg_qual = (float(sum(a)) / float(total_clusters)) * 100
        return avg_qual

    def get_qual_graph(self, qualitymetrics, avg_qual):
        quality_df = qualitymetrics.df

        cols_to_drop = ['cycle', 'lane', 'tile', 'total']
        quality_df = quality_df.drop(cols_to_drop, axis=1)
        qual_scores = quality_df.sum(axis=0).values
        qual_list = map(int, qual_scores.tolist())
        less_than_q30 = qual_list[0: 29]
        more_than_q30 = qual_list[30: 50]

        myint = 1000000
        y1 = [x / myint for x in less_than_q30]
        y2 = [x / myint for x in more_than_q30]
        x1 = range(1, 30)
        x2 = range(30, 50)
        y_max = max(y1 + y2)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.bar(x1, y1, color='b', width=1.0)
        ax.bar(x2, y2, color='g', width=1.0)
        ax.set_ylabel('Total (millions)')
        ax.set_xlabel('Q score')
        ax.plot([30, 30], [0, y_max], color='g', linewidth=1, zorder=5)
        percent = '%s%%' % format(avg_qual, '.2f')
        ax.text(23, y_max - 100, percent, fontsize=16, color='g')

    def get_qual_heatmap(self, qualitymetrics):
        quality_df = qualitymetrics.df
        cols_to_drop = ['lane', 'tile', 'total']
        quality_df = quality_df.drop(cols_to_drop, axis=1)
        quality_df.columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                              26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                              48, 49, 50, 'cycle']
        max_cycle = quality_df['cycle'].max()
        grouped_by_cycle = quality_df.groupby(['cycle'])
        x = grouped_by_cycle.aggregate(np.sum)
        x = x.transpose()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        my_cmap = plt.cm.get_cmap('GnBu')
        my_cmap.set_under(color='white')
        ax.pcolor(x, cmap=my_cmap, vmin=0.0001)

        major_xticks = np.arange(0, max_cycle + 20, 20)
        major_yticks = np.arange(0, 50, 10)
        ax.set_xticks(major_xticks)
        ax.set_yticks(major_yticks)
        ax.grid(b=True, which='both', color='0.85', linestyle='-')
        ax.spines['right'].set_color('0.85')
        ax.spines['top'].set_color('0.85')
        ax.spines['bottom'].set_color('0.85')
        ax.set_xlabel('Cycle')
        ax.set_ylabel('Q Score')

    def get_clusters_per_sample(self, indexmetrics, tilemetrics):
        total_aligned_clusters = float(tilemetrics.num_clusters * tilemetrics.aligned)
        pf_aligned_clusters = float(indexmetrics.df['clusters'].sum())
        percent_pf = format(float(pf_aligned_clusters / total_aligned_clusters) * 100, '.4f')

        group_by_index = indexmetrics.df.groupby(['name_str', 'index_str'], sort=False, as_index=False)
        grouped_df = group_by_index.aggregate(np.sum)
        cols_to_drop = ['lane', 'read', 'tile']
        grouped_df = grouped_df.drop(cols_to_drop, axis=1)
        sample_id = grouped_df['name_str'].tolist()
        clusters = grouped_df['clusters'].tolist()
        s = grouped_df['index_str'].str.split('-', expand=True)
        index1 = s[0]
        index2 = s[1]

        percent_clusters = [format((float(x) / total_aligned_clusters) * 100, '.4f') for x in clusters]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(sample_id, percent_clusters)
        ax.set_ylabel('% Reads')
        ax.set_xlabel('Sample_Id')

        return sample_id, index1, index2, percent_clusters, percent_pf, total_aligned_clusters, pf_aligned_clusters
