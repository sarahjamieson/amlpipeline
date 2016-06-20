from pylatex import Document, Section, Tabular, Package, Command, Figure, SubFigure
from pylatex.utils import NoEscape, bold
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np


class CreatePDF(object):
    def __init__(self, tilemetrics, controlmetrics, errormetrics, extractionmetrics, indexmetrics, qualitymetrics,
                 corintmetrics, worksheet):
        self.tile = tilemetrics
        self.control = controlmetrics
        self.error = errormetrics
        self.extract = extractionmetrics
        self.index = indexmetrics
        self.quality = qualitymetrics
        self.corint = corintmetrics
        self.worksheet = worksheet

    def create_pdf(self):
        """Creates a LaTeX PDF using the Python module PyLatex (https://jeltef.github.io/PyLaTeX/latest/).

            Note: LaTeX must be installed in /usr/local/. Development employed TexLive 2015.


        """
        pdflatex = '/usr/local/texlive/2015/bin/x86_64-linux/pdflatex'
        doc = Document()
        doc.packages.append(Package('geometry', options=['tmargin=0.75in', 'lmargin=0.75in', 'rmargin=0.75in']))
        doc.packages.append(Package('datetime', options=['ddmmyyyy']))
        doc.packages.append(Package('needspace'))
        doc.preamble.append(Command('newdateformat',
                                    NoEscape(r'mydate}{\twodigit{\THEDAY}/\twodigit{\THEMONTH}/\THEYEAR')))
        doc.append(Command('begin', 'center'))
        doc.append(Command('Large', bold('MiSeq Quality Check')))
        doc.append(Command('end', 'center'))
        doc.append(Command('begin', 'flushright'))
        doc.append(Command('Large', '<Run_name>'))
        doc.append(Command('end', 'flushright'))
        doc.append(Command('begin', 'flushright'))
        doc.append(Command('Large', NoEscape(r'\mydate\today')))
        doc.append(Command('end', 'flushright'))

        avg_qual = self.get_avg_qual()
        self.get_qual_graph(avg_qual)
        first_read, second_read = self.get_avg_qual_per_read()
        doc.append(Command('needspace', '20em'))
        with doc.create(Section('Quality data')):
            with doc.create(Tabular(NoEscape(r'p{5cm}|c|c'))) as table:
                table.add_row(('Data', 'Value', 'Pass/Fail'))
                table.add_hline()
                table.add_row(
                    ('Mean Cluster Density (k/mm2)', format(self.tile.mean_cluster_density / 1000, '.2f'), ''))
                table.add_row(('Clusters passed filter (%)', '%s' % (format(self.tile.percent_pf_clusters, '.2f')),
                               ''))
                table.add_row(('Average >= Q30', '%s' % (format(avg_qual, '.2f')), ''))
                table.add_row(('1st full read >= Q30', '%s' % (format(first_read, '.2f')), ''))
                table.add_row(('2nd full read >= Q30', '%s' % (format(second_read, '.2f')), ''))

            with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))):
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_plot()
                    plot.add_caption('Q-score distribution plot (all reads all cycles)')
                self.get_qual_heatmap()
                with doc.create(SubFigure()) as plot:
                    plot.add_plot()
                    plot.add_caption('Q-score heat map')
                self.get_clusters_heatmap()
                with doc.create(SubFigure()) as plot:
                    plot.add_plot()
                    plot.add_caption('Cluster density per tile')
        read_1_phas, read_1_prephas, read_2_phas, read_2_prephas = self.get_phas_prephas()
        with doc.create(Section('Phas/Prephas data')):
            with doc.create(Tabular(NoEscape(r'p{5cm}|c|c'))) as table:
                table.add_row(('Data', 'Value', 'Pass/Fail'))
                table.add_hline()
                table.add_row(
                    ('1st full read', '%s / %s' % (format(read_1_phas, '.3f'), format(read_1_prephas, '.3f')), ''))
                table.add_row(
                    ('2nd full read', '%s / %s' % (format(read_2_phas, '.3f'), format(read_2_prephas, '.3f')), ''))
        sample_id, index1, index2, percent_clusters, percent_pf, total_aligned_clusters, pf_aligned_clusters = \
            self.get_indexing()
        total_samples = len(sample_id)

        doc.append(Command('needspace', '10em'))
        with doc.create(Section('Indexing')):
            doc.append(Command('begin', 'center'))
            with doc.create(Tabular(NoEscape(r'c|c|c|c|c'))) as table:
                table.add_row(('Total Reads', 'PF Reads', '% Reads Identified (PF)', 'Min', 'Max'))
                table.add_hline()
                table.add_row(
                    ('%s' % int(total_aligned_clusters), '%s' % int(pf_aligned_clusters), '%s' % percent_pf,
                     '%s' % min(percent_clusters), '%s' % max(percent_clusters)))
            doc.append(Command('end', 'center'))
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

        doc.generate_pdf('%s_InterOp_Results' % self.worksheet, clean_tex=False, compiler=pdflatex)
        # os.system("xdg-open /home/cuser/PycharmProjects/amlpipeline/output.pdf")

    def get_avg_qual(self):
        """Calculates the percentage of clusters with an average quality score of >=30 using InterOp qualitymetrics.

        :return: avg_qual: percentage of clusters.
        """
        # 1) Extract data into pandas data frame
        quality_df = self.quality.df

        # 2) Calculate number of clusters in each row, add to new column named "Total" and sum to give overall total.
        col_list = range(0, 50)
        qual_data = quality_df[col_list].sum(axis=1)
        quality_df['total'] = qual_data.tolist()
        total_clusters = quality_df['total'].sum(axis=0)

        # 3) Calculate number of clusters with q>=30 and divide by total clusters to generate %.
        col_list = range(29, 50)
        qual_over_30 = quality_df[col_list].sum(axis=1)
        avg_qual = (float(sum(qual_over_30)) / float(total_clusters)) * 100

        return avg_qual

    def get_avg_qual_per_read(self):
        """Calculates percentage of clusters with average quality score of >=30 per read using InterOp qualitymetrics.

        :return: first_read_avg: percentage for first full read (~9-159 cycles)
        :return: second_read_avg: percentage for second full read (~168-318 cycles)
        """
        # 1) Extract data into pandas data frame
        quality_df = self.quality.df

        # 2) Group quality data by cycle
        grouped_by_cycle = quality_df.groupby(['cycle'])
        quality_by_cycle_df = grouped_by_cycle.aggregate(np.sum)

        # 3) Calculate number of cycles
        max_cycle = quality_df['cycle'].max()

        # 4) Calculate number of clusters per cycle and then per read
        total_clusters_per_cycle = quality_by_cycle_df.iloc[0]['total']
        clusters_per_read = total_clusters_per_cycle * (max_cycle - 8)

        # 5) Calculate number of clusters with q score >= 30
        col_list = range(29, 50)
        q_over_30 = quality_by_cycle_df[col_list].sum(axis=1)
        q_over_30_list = q_over_30.tolist()

        # 6) Split list to separate cycles in first read and second read.
        first_read_list = q_over_30_list[8:max_cycle / 2]
        second_read_list = q_over_30_list[(max_cycle / 2) + 8:max_cycle + 1]

        # 7) Sum of all clusters in each read and divide by total to generate %.
        first_read_avg = (float(sum(first_read_list)) / float(clusters_per_read)) * 100
        second_read_avg = (float(sum(second_read_list)) / float(clusters_per_read)) * 100

        return first_read_avg, second_read_avg

    def calculate_percentage(self, value):
        """Takes a value and divides it by the total number of clusters passing filter (PF) to generate %.

        :param value
        :return: percentage
        """
        return (value / self.tile.num_clusters_pf) * 100

    def get_qual_graph(self, avg_qual):
        """Takes

        :param avg_qual:
        :return:
        """
        quality_df = self.quality.df

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

    def get_qual_heatmap(self):
        quality_df = self.quality.df
        cols_to_drop = ['lane', 'tile', 'total']
        quality_df = quality_df.drop(cols_to_drop, axis=1)
        quality_df.columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                              26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                              48, 49, 50, 'cycle']
        max_cycle = quality_df['cycle'].max()
        grouped_by_cycle = quality_df.groupby(['cycle'])
        x = grouped_by_cycle.aggregate(np.sum)
        x = x.transpose()
        x = x.apply(self.calculate_percentage)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        colours = [[255, 255, 255], [252, 255, 252], [248, 253, 248], [244, 252, 244], [240, 251, 240], [236, 250, 236],
                   [232, 248, 232], [228, 247, 228], [224, 246, 224], [220, 245, 220], [216, 244, 216], [212, 242, 212],
                   [208, 241, 208], [204, 240, 204], [200, 239, 200], [196, 237, 196], [192, 236, 192], [188, 235, 188],
                   [184, 234, 184], [180, 233, 180], [176, 231, 176], [172, 230, 172], [168, 229, 168], [163, 228, 163],
                   [159, 226, 159], [155, 225, 155], [151, 224, 151], [147, 223, 147], [143, 221, 143], [139, 220, 139],
                   [135, 219, 135], [131, 218, 131], [127, 217, 127], [123, 215, 123], [119, 214, 119], [115, 213, 115],
                   [111, 212, 111], [107, 210, 107], [103, 209, 103], [99, 208, 99], [95, 207, 95], [91, 206, 91],
                   [87, 204, 87], [83, 203, 83], [79, 202, 79], [75, 201, 75], [71, 199, 71], [67, 198, 67],
                   [63, 197, 63], [59, 196, 59], [55, 195, 55], [51, 193, 51], [47, 192, 47], [43, 191, 43],
                   [39, 190, 39], [35, 188, 35], [31, 187, 31], [27, 186, 27], [23, 185, 23], [19, 184, 19],
                   [15, 182, 15], [11, 181, 11], [7, 180, 7], [3, 179, 3], [1, 178, 0], [5, 179, 0], [9, 180, 0],
                   [13, 182, 0], [17, 183, 0], [21, 184, 0], [25, 185, 0], [29, 186, 0], [33, 188, 0], [37, 189, 0],
                   [41, 190, 0], [45, 191, 0], [49, 192, 0], [53, 194, 0], [57, 195, 0], [61, 196, 0], [65, 197, 0],
                   [69, 198, 0], [73, 200, 0], [77, 201, 0], [81, 202, 0], [85, 203, 0], [89, 204, 0], [93, 206, 0],
                   [97, 207, 0], [101, 208, 0], [105, 209, 0], [109, 210, 0], [113, 212, 0], [117, 213, 0],
                   [121, 214, 0], [125, 215, 0], [129, 216, 0], [133, 218, 0], [137, 219, 0], [141, 220, 0],
                   [145, 221, 0], [149, 222, 0], [153, 224, 0], [157, 225, 0], [161, 226, 0], [165, 227, 0],
                   [169, 228, 0], [173, 230, 0], [177, 231, 0], [181, 232, 0], [185, 233, 0], [189, 234, 0],
                   [193, 236, 0], [197, 237, 0], [201, 238, 0], [205, 239, 0], [209, 240, 0], [213, 242, 0],
                   [217, 243, 0], [221, 244, 0], [225, 245, 0], [229, 246, 0], [233, 248, 0], [237, 249, 0],
                   [241, 250, 0], [245, 251, 0], [249, 252, 0], [253, 254, 0], [255, 254, 0], [255, 252, 1],
                   [255, 251, 1], [255, 250, 1], [255, 249, 2], [254, 247, 2], [254, 246, 3], [254, 245, 3],
                   [254, 244, 4], [254, 243, 4], [254, 241, 4], [254, 240, 5], [254, 239, 5], [254, 238, 6],
                   [253, 236, 6], [253, 235, 7], [253, 234, 7], [253, 233, 7], [253, 232, 8],
                   [253, 230, 8], [253, 229, 9], [253, 228, 9], [253, 227, 10], [252, 225, 10], [252, 224, 10],
                   [252, 223, 11], [252, 222, 11], [252, 221, 12], [252, 219, 12], [252, 218, 13], [252, 217, 13],
                   [252, 216, 13], [251, 214, 14], [251, 213, 14], [251, 212, 15], [251, 211, 15], [251, 210, 15],
                   [251, 208, 16], [251, 207, 16], [251, 206, 17], [251, 205, 17], [250, 203, 18], [250, 202, 18],
                   [250, 201, 18], [250, 200, 19], [250, 199, 19], [250, 197, 20], [250, 196, 20], [250, 195, 21],
                   [250, 194, 21], [249, 192, 21], [249, 191, 22], [249, 190, 22], [249, 189, 23], [249, 188, 23],
                   [249, 186, 24], [249, 185, 24], [249, 184, 24], [249, 183, 25], [249, 181, 25], [248, 180, 26],
                   [248, 179, 26], [248, 178, 27], [248, 176, 27], [248, 174, 27], [248, 171, 26], [248, 169, 26],
                   [249, 166, 25], [249, 163, 25], [249, 160, 25], [249, 158, 24], [249, 155, 24], [249, 152, 23],
                   [249, 149, 23], [249, 146, 22], [249, 144, 22], [250, 141, 22], [250, 138, 21], [250, 135, 21],
                   [250, 133, 20], [250, 130, 20], [250, 127, 20], [250, 124, 19], [250, 122, 19], [251, 119, 18],
                   [251, 116, 18], [251, 113, 17], [251, 111, 17], [251, 108, 17], [251, 105, 16], [251, 102, 16],
                   [251, 99, 15], [252, 97, 15], [252, 94, 14], [252, 91, 14], [252, 88, 14], [252, 86, 13],
                   [252, 83, 13], [252, 80, 12], [252, 77, 12], [253, 75, 11], [253, 72, 11], [253, 69, 11],
                   [253, 66, 10], [253, 64, 10], [253, 61, 9], [253, 58, 9], [253, 55, 8], [253, 53, 8], [254, 50, 8],
                   [254, 47, 7], [254, 44, 7], [254, 41, 6], [254, 39, 6], [254, 36, 6], [254, 33, 5], [254, 30, 5],
                   [255, 28, 4], [255, 25, 4], [255, 22, 3], [255, 19, 3], [255, 17, 3], [255, 14, 2], [255, 11, 2],
                   [255, 8, 1], [255, 6, 1], [255, 3, 0], [255, 0, 0]]

        colour_codes = [[float(j) / 255 for j in i] for i in colours]
        my_cmap = ListedColormap(colour_codes)
        heatmap = ax.pcolor(x, cmap=my_cmap)

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
        plt.colorbar(heatmap)

    def get_phas_prephas(self):
        read_1_phas = self.tile.mean_phasing[0] * 100
        read_2_phas = self.tile.mean_phasing[1] * 100
        read_1_prephas = self.tile.mean_prephasing[0] * 100
        read_2_prephas = self.tile.mean_prephasing[1] * 100

        return read_1_phas, read_1_prephas, read_2_phas, read_2_prephas

    def get_indexing(self):
        total_aligned_clusters = float(self.tile.num_clusters * self.tile.aligned)
        pf_aligned_clusters = float(self.index.df['clusters'].sum())
        percent_pf = format(float(pf_aligned_clusters / total_aligned_clusters) * 100, '.4f')

        group_by_index = self.index.df.groupby(['name_str', 'index_str'], sort=False, as_index=False)
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

    def get_clusters_heatmap(self):
        tile_df = self.tile.df[self.tile.df['code'] == 100]
        cols_to_drop = ['code', 'lane']
        tile_df = tile_df.drop(cols_to_drop, axis=1)
        tile_df = tile_df.sort_values('tile')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        my_cmap = plt.cm.get_cmap('jet')

        num_tiles = self.tile.num_tiles
        tiles_11 = tile_df['value'][0:num_tiles / 2].tolist()
        tiles_21 = tile_df['value'][num_tiles / 2:num_tiles].tolist()
        myint = 1000
        tiles_11_div = [x / myint for x in tiles_11]
        tiles_21_div = [x / myint for x in tiles_21]

        data = []
        item = 0
        while item < (num_tiles / 2):
            data.append([tiles_11_div[item], tiles_21_div[item]])
            item += 1

        heatmap = ax.pcolor(data, cmap=my_cmap, vmin=0, vmax=2400)

        x_ticks = [1, 2]
        x_labels = ['11', '21']
        y_labels = []
        y_ticks = []
        y = 1
        while y <= (num_tiles / 2):
            if y < 10:
                y_labels.append('0%s' % str(y))
            else:
                y_labels.append(str(y))
            y_ticks.append(float(y))
            y += 1

        y_ticks_cen = [x - 0.5 for x in y_ticks]
        x_ticks_cen = [x - 0.5 for x in x_ticks]

        ax.set_xticks(x_ticks_cen, minor=False)
        ax.set_xticklabels(x_labels, ha='center')

        ax.set_yticks(y_ticks_cen, minor=False)
        ax.set_yticklabels(y_labels, minor=False)

        plt.colorbar(heatmap)
