import pandas as pd
from zipfile import ZipFile
from pylatex import Document, Section, Tabular, Package, Command, Figure, SubFigure, Description
from pylatex.utils import NoEscape, bold
import os


class CreateFastQCPDF(object):
    def __init__(self, fastqc):
        self.fastqc = fastqc

    def create_pdf(self):
        # 1) Extract files from zip folder
        with ZipFile(self.fastqc, 'w') as myzip:
            myzip.extractall()

        # 2) Put summary.txt and fastqc_data.txt files into dataframes
        summary_df = pd.read_table('04-R1.qfilter_fastqc/summary.txt', header=None, names=['Score', 'Parameter'],
                                   usecols=[0, 1])
        score_list = summary_df['Score'].tolist()  # not currently used, may be needed
        parameter_list = summary_df['Parameter'].tolist()
        summary_dict = dict(zip(parameter_list, score_list))

        basic_stats_df = pd.read_table('04-R1.qfilter_fastqc/fastqc_data.txt', header=None, names=['Property', 'Value'],
                                       usecols=[0, 1], skiprows=3, nrows=7)

        # 3) Create PDF with basic statistics, summary data and png images
        doc = Document()
        doc.packages.append(Package('geometry', options=['tmargin=0.75in', 'lmargin=0.75in', 'rmargin=0.75in']))
        doc.packages.append(Package('subcaption'))
        doc.packages.append(Package('xcolor'))
        doc.append(Command(NoEscape(r'renewcommand{\baselinestretch}'), '1.0'))
        doc.append(Command('begin', 'center'))
        doc.append(Command('Large', bold('FastQC Quality Results')))
        doc.append(Command('end', 'center'))

        with doc.create(Section('Basic Statistics')):
            with doc.create(Description()) as desc:
                for row_index, row in basic_stats_df.iterrows():
                    desc.add_item("%s:" % row['Property'], "%s" % row['Value'])

        with doc.create(Section('Summary')):
            with doc.create(Tabular(NoEscape(r'p{5cm}|c'))) as table:
                table.add_row(('Parameter', 'Pass/Warning/Fail'))
                table.add_hline()
                for row_index, row in summary_df.iterrows():
                    if row['Score'] == 'PASS':
                        colour = 'green'
                    elif row['Score'] == 'WARN':
                        colour = 'orange'
                    elif row['Score'] == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    table.add_row(('%s' % row['Parameter'], NoEscape(r'\textcolor{%s}{%s}' % (colour, row['Score']))))

            with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))):
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/per_base_quality.png')
                    plot.add_caption('Per base sequence quality')
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/per_tile_quality.png')
                    plot.add_caption('Per tile sequence quality')
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/per_sequence_quality.png')
                    plot.add_caption('Per sequence quality scores')
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/per_base_sequence_content.png')
                    plot.add_caption('Per base sequence content')
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/per_sequence_gc_content.png')
                    plot.add_caption('Per sequence GC content')
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/per_base_n_content.png')
                    plot.add_caption('Per base N content')
            with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))):
                doc.append(Command('ContinuedFloat'))
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/sequence_length_distribution.png')
                    plot.add_caption('Sequence Length Distribution')
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/duplication_levels.png')
                    plot.add_caption('Sequence Duplication Levels')
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/adapter_content.png')
                    plot.add_caption('Adapter content: %s')
                with doc.create(SubFigure()) as plot:
                    plot.add_image('04-R1.qfilter_fastqc/Images/kmer_profiles.png')
                    plot.add_caption('Kmer content: %s')

        pdflatex = '/usr/local/texlive/2015/bin/x86_64-linux/pdflatex'
        doc.generate_pdf('fastqc', clean_tex=False, compiler=pdflatex)
        os.system('mv /home/cuser/PycharmProjects/amlpipeline/fastqc.pdf /media/sf_sarah_share/MiSeq_quality_outputs/')
