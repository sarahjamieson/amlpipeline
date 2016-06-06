import pandas as pd
from zipfile import ZipFile
from pylatex import Document, Section, Tabular, Package, Command, Figure, SubFigure, Description
from pylatex.utils import NoEscape, bold
import os


class CreateFastQCPDF(object):
    def __init__(self, sample):
        self.sample = sample

    def create_pdf(self):
        # R1
        # ----------------------------------------------------------------------------------
        # 2) Put summary.txt and fastqc_data.txt files into dataframes then dictionaries
        # i) Trimmed data
        summary_df_trim = pd.read_table('%sR1_001.qfilter_fastqc/summary.txt' % self.sample, header=None,
                                        names=['Score', 'Parameter'], usecols=[0, 1])
        score_list_trim = summary_df_trim['Score'].tolist()  # not currently used, may be needed
        parameter_list_trim = summary_df_trim['Parameter'].tolist()
        sum_dict_trim = dict(zip(parameter_list_trim, score_list_trim))

        # ii) Original non-trimmed data
        summary_df = pd.read_table('%sR1_001_fastqc/summary.txt' % self.sample, header=None,
                                   names=['Score', 'Parameter'], usecols=[0, 1])
        score_list = summary_df_trim['Score'].tolist()  # not currently used, may be needed
        parameter_list = summary_df_trim['Parameter'].tolist()
        sum_dict = dict(zip(parameter_list, score_list))

        # basic stats files
        # Trimmed data
        basic_stats_df_trim = pd.read_table('%sR1_001.qfilter_fastqc/fastqc_data.txt' % self.sample, header=None,
                                            names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
        property_list_trim = basic_stats_df_trim['Property'].tolist()
        value_list_trim = basic_stats_df_trim['Value'].tolist()
        stats_dict_trim = dict(zip(property_list_trim, value_list_trim))
        # Original, non-trimmed data
        basic_stats_df = pd.read_table('%sR1_001_fastqc/fastqc_data.txt' % self.sample, header=None,
                                       names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
        property_list = basic_stats_df['Property'].tolist()
        value_list = basic_stats_df['Value'].tolist()
        stats_dict = dict(zip(property_list, value_list))

        # 3) Create PDF with basic statistics, summary data and png images
        doc = Document()
        doc.packages.append(Package('geometry', options=['tmargin=0.75in', 'lmargin=0.5in', 'rmargin=0.5in']))
        doc.packages.append(Package('subcaption'))
        doc.packages.append(Package('xcolor'))
        doc.append(Command('makeatletter'))
        doc.append(Command('setlength', NoEscape(r'\@fptop}{0pt')))
        doc.append(Command('makeatother'))
        doc.append(Command(NoEscape(r'renewcommand{\baselinestretch}'), '1.0'))
        doc.append(Command('begin', 'center'))
        doc.append(Command('Large', bold('Sample Quality Results')))
        doc.append(Command('end', 'center'))

        with doc.create(Section('Basic Statistics')):
            with doc.create(Description()) as desc:
                desc.add_item("Filename:", "%s" % stats_dict.get('Filename'))
                desc.add_item("File type:", "%s" % stats_dict.get('File type'))
                desc.add_item("Encoding:", "%s" % stats_dict.get('Encoding'))
            with doc.create(Tabular(NoEscape(r'p{5.5cm}|c|c'))) as table:
                table.add_row(('', 'Before trimming', 'After trimming'))
                table.add_hline()
                table.add_row(('Total Sequences', '%s' % stats_dict.get('Total Sequences'),
                               '%s' % stats_dict_trim.get('Total Sequences')))
                table.add_row(('Sequences flagged as poor quality',
                               '%s' % stats_dict.get('Sequences flagged as poor quality'),
                               '%s' % stats_dict_trim.get('Sequences flagged as poor quality')))
                table.add_row(('Sequence length', '%s' % stats_dict.get('Sequence length'),
                               '%s' % stats_dict_trim.get('Sequence length')))
                table.add_row(('%GC', '%s' % stats_dict.get('%GC'), '%s' % stats_dict_trim.get('%GC')))

        with doc.create(Section('FastQC')):
            with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))):
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%s/Images/per_base_quality.png' % self.fastqc)
                    plot.add_caption('Per base sequence quality BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    if sum_dict_trim.get('Per base sequence quality') == 'PASS':
                        colour = 'green'
                    elif sum_dict_trim.get('Per base sequence quality') == 'WARN':
                        colour = 'orange'
                    elif sum_dict_trim.get('Per base sequence quality') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%s/Images/per_base_quality.png' % self.fastqc_trim)
                    plot.add_caption(NoEscape(r'Per base sequence quality AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, sum_dict_trim.get('Per base sequence quality'))))
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '5 mm'))
                    plot.add_image('%s/Images/per_sequence_gc_content.png' % self.fastqc)
                    plot.add_caption('Per sequence GC content BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '5 mm'))
                    if sum_dict_trim.get('Per sequence GC content') == 'PASS':
                        colour = 'green'
                    elif sum_dict_trim.get('Per sequence GC content') == 'WARN':
                        colour = 'orange'
                    elif sum_dict_trim.get('Per sequence GC content') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%s/Images/per_sequence_gc_content.png' % self.fastqc_trim)
                    plot.add_caption(NoEscape(r'Per sequence GC content AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, sum_dict_trim.get('Per sequence GC content'))))
            with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))):
                doc.append(Command('ContinuedFloat'))
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%s/Images/sequence_length_distribution.png' % self.fastqc)
                    plot.add_caption('Sequence Length Distribution BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('hspace', '10 mm'))
                    if sum_dict_trim.get('Sequence Length Distribution') == 'PASS':
                        colour = 'green'
                    elif sum_dict_trim.get('Sequence Length Distribution') == 'WARN':
                        colour = 'orange'
                    elif sum_dict_trim.get('Sequence Length Distribution') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%s/Images/sequence_length_distribution.png' % self.fastqc_trim)
                    plot.add_caption(NoEscape(r'Sequence Length Distribution AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, sum_dict_trim.get('Sequence Length Distribution'))))
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '10 mm'))
                    plot.add_image('%s/Images/adapter_content.png' % self.fastqc)
                    plot.add_caption('Adapter content BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '10 mm'))
                    doc.append(Command('hspace', '10 mm'))
                    if sum_dict_trim.get('Adapter Content') == 'PASS':
                        colour = 'green'
                    elif sum_dict_trim.get('Adapter Content') == 'WARN':
                        colour = 'orange'
                    elif sum_dict_trim.get('Adapter Content') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%s/Images/adapter_content.png' % self.fastqc_trim)
                    plot.add_caption(NoEscape(r'Adapter content AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, sum_dict_trim.get('Adapter Content'))))

        with doc.create(Section('BamStats')):
            with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))):
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%s.bwa.drm.sorted.bam.stats-quals-hm.png' % self.name)
                    plot.add_caption('Base quality per cycle')
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('hspace', '10 mm'))
                    plot.add_image('%s.bwa.drm.sorted.bam.stats-insert-size.png' % self.name)
                    plot.add_caption('Fragment size')
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '10 mm'))
                    plot.add_image('%s.bwa.drm.sorted.bam.stats-quals2.png' % self.name)
                    plot.add_caption('Quality per cycle')

        pdflatex = '/usr/local/texlive/2015/bin/x86_64-linux/pdflatex'
        doc.generate_pdf('%s' % self.name, clean_tex=False, compiler=pdflatex)
        os.system('mv /home/cuser/PycharmProjects/amlpipeline/%s.pdf /media/sf_sarah_share/MiSeq_quality_outputs/'
                  % self.name)

pdf = CreateFastQCPDF('04-D15-22373-HT-Nextera-Myeloid-Val1-Repeat_S4_L001_R1_001_fastqc',
                      '04-D15-22373-HT-Nextera-Myeloid-Val1-Repeat_S4_L001_R1_001.qfilter_fastqc')
pdf.create_pdf()
