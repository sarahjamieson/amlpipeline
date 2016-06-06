import pandas as pd
from pylatex import Document, Section, Tabular, Package, Command, Figure, SubFigure, Description, Subsection
from pylatex.utils import NoEscape, bold
import os


class CreateFastQCPDF(object):
    def __init__(self, sample):
        self.sample = sample

    def create_pdf(self):
        r1_summary_trim_dict, r1_stats_trim_dict, r2_summary_trim_dict, r2_stats_trim_dict = self.get_trimmed_data()
        r1_summary_dict, r1_stats_dict, r2_summary_dict, r2_stats_dict = self.get_original_data()

        doc = Document()
        doc.packages.append(Package('geometry', options=['margin=0.75in']))
        doc.packages.append(Package('subcaption'))
        doc.packages.append(Package('xcolor'))
        doc.packages.append(Package('needspace'))
        doc.append(Command('makeatletter'))
        doc.append(Command('setlength', NoEscape(r'\@fptop}{0pt')))
        doc.append(Command('makeatother'))
        doc.append(Command(NoEscape(r'renewcommand{\baselinestretch}'), '1.0'))
        doc.append(Command('begin', 'center'))
        doc.append(Command('Large', bold('Sample Quality Results')))
        doc.append(Command('end', 'center'))

        with doc.create(Section('Basic Statistics')):
            with doc.create(Description()) as desc:
                desc.add_item("Sample:", "%s" % r1_stats_dict.get('Filename')[:-16])
                desc.add_item("File type:", "%s" % r1_stats_dict.get('File type'))
                desc.add_item("Encoding:", "%s" % r1_stats_dict.get('Encoding'))
            with doc.create(Tabular(NoEscape(r'p{5.5cm}|c|c|c|c'))) as table:
                table.add_row(('', 'R1', '', 'R2', ''))
                table.add_hline()
                table.add_row(('', 'Before trimming', 'After trimming', 'Before trimming', 'After trimming'))
                table.add_hline()
                table.add_row(('Total Sequences',
                               '%s' % r1_stats_dict.get('Total Sequences'),
                               '%s' % r1_stats_trim_dict.get('Total Sequences'),
                               '%s' % r2_stats_dict.get('Total Sequences'),
                               '%s' % r2_stats_trim_dict.get('Total Sequences')))
                table.add_row(('Sequences flagged as poor quality',
                               '%s' % r1_stats_dict.get('Sequences flagged as poor quality'),
                               '%s' % r1_stats_trim_dict.get('Sequences flagged as poor quality'),
                               '%s' % r2_stats_dict.get('Sequences flagged as poor quality'),
                               '%s' % r2_stats_trim_dict.get('Sequences flagged as poor quality')))
                table.add_row(('Sequence length',
                               '%s' % r1_stats_dict.get('Sequence length'),
                               '%s' % r1_stats_trim_dict.get('Sequence length'),
                               '%s' % r2_stats_dict.get('Sequence length'),
                               '%s' % r2_stats_trim_dict.get('Sequence length')))
                table.add_row(('%GC',
                               '%s' % r1_stats_dict.get('%GC'),
                               '%s' % r1_stats_trim_dict.get('%GC'),
                               '%s' % r2_stats_dict.get('%GC'),
                               '%s' % r2_stats_trim_dict.get('%GC')))

        with doc.create(Section('FastQC')):
            with doc.create(Figure(position='!htb', placement=NoEscape(r'\centering'))) as fig:
                fig.add_caption('Per base sequence quality')
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%sR1_001_fastqc/Images/per_base_quality.png' % self.sample)
                    plot.add_caption('R1 BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    if r1_summary_trim_dict.get('Per base sequence quality') == 'PASS':
                        colour = 'green'
                    elif r1_summary_trim_dict.get('Per base sequence quality') == 'WARN':
                        colour = 'orange'
                    elif r1_summary_trim_dict.get('Per base sequence quality') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%sR1_001.qfilter_fastqc/Images/per_base_quality.png' % self.sample)
                    plot.add_caption(NoEscape(r'R1 AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, r1_summary_trim_dict.get('Per base sequence quality'))))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%sR2_001_fastqc/Images/per_base_quality.png' % self.sample)
                    plot.add_caption('R2 BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    if r2_summary_trim_dict.get('Per base sequence quality') == 'PASS':
                        colour = 'green'
                    elif r2_summary_trim_dict.get('Per base sequence quality') == 'WARN':
                        colour = 'orange'
                    elif r2_summary_trim_dict.get('Per base sequence quality') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%sR2_001.qfilter_fastqc/Images/per_base_quality.png' % self.sample)
                    plot.add_caption(NoEscape(r'R2 AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, r2_summary_trim_dict.get('Per base sequence quality'))))

            with doc.create(Figure(position='!htb', placement=NoEscape(r'\centering'))) as fig:
                fig.add_caption('Per sequence GC content')
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '5 mm'))
                    plot.add_image('%sR1_001_fastqc/Images/per_sequence_gc_content.png' % self.sample)
                    plot.add_caption('R1 BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '5 mm'))
                    if r1_summary_trim_dict.get('Per sequence GC content') == 'PASS':
                        colour = 'green'
                    elif r1_summary_trim_dict.get('Per sequence GC content') == 'WARN':
                        colour = 'orange'
                    elif r1_summary_trim_dict.get('Per sequence GC content') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%sR1_001.qfilter_fastqc/Images/per_sequence_gc_content.png' % self.sample)
                    plot.add_caption(NoEscape(r'R1 AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, r1_summary_trim_dict.get('Per sequence GC content'))))
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '5 mm'))
                    plot.add_image('%sR2_001_fastqc/Images/per_sequence_gc_content.png' % self.sample)
                    plot.add_caption('R2 BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '5 mm'))
                    if r2_summary_trim_dict.get('Per sequence GC content') == 'PASS':
                        colour = 'green'
                    elif r2_summary_trim_dict.get('Per sequence GC content') == 'WARN':
                        colour = 'orange'
                    elif r2_summary_trim_dict.get('Per sequence GC content') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%sR2_001.qfilter_fastqc/Images/per_sequence_gc_content.png' % self.sample)
                    plot.add_caption(NoEscape(r'R2 AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, r2_summary_trim_dict.get('Per sequence GC content'))))

            with doc.create(Figure(position='!htb', placement=NoEscape(r'\centering'))) as fig:
                fig.add_caption('Sequence Length Distribution')
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%sR1_001_fastqc/Images/sequence_length_distribution.png' % self.sample)
                    plot.add_caption('R1 BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    if r1_summary_trim_dict.get('Sequence Length Distribution') == 'PASS':
                        colour = 'green'
                    elif r1_summary_trim_dict.get('Sequence Length Distribution') == 'WARN':
                        colour = 'orange'
                    elif r1_summary_trim_dict.get('Sequence Length Distribution') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%sR1_001.qfilter_fastqc/Images/sequence_length_distribution.png' % self.sample)
                    plot.add_caption(NoEscape(r'R1 AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, r1_summary_trim_dict.get('Sequence Length Distribution'))))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%sR2_001_fastqc/Images/sequence_length_distribution.png' % self.sample)
                    plot.add_caption('R2 BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    if r2_summary_trim_dict.get('Sequence Length Distribution') == 'PASS':
                        colour = 'green'
                    elif r2_summary_trim_dict.get('Sequence Length Distribution') == 'WARN':
                        colour = 'orange'
                    elif r2_summary_trim_dict.get('Sequence Length Distribution') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%sR2_001.qfilter_fastqc/Images/sequence_length_distribution.png' % self.sample)
                    plot.add_caption(NoEscape(r'R2 AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, r2_summary_trim_dict.get('Sequence Length Distribution'))))

            with doc.create(Figure(position='!htb', placement=NoEscape(r'\centering'))) as fig:
                fig.add_caption('Adapter Content')
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%sR1_001_fastqc/Images/adapter_content.png' % self.sample)
                    plot.add_caption('R1 BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    if r1_summary_trim_dict.get('Adapter Content') == 'PASS':
                        colour = 'green'
                    elif r1_summary_trim_dict.get('Adapter Content') == 'WARN':
                        colour = 'orange'
                    elif r1_summary_trim_dict.get('Adapter Content') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%sR1_001.qfilter_fastqc/Images/adapter_content.png' % self.sample)
                    plot.add_caption(NoEscape(r'R1 AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, r1_summary_trim_dict.get('Adapter Content'))))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%sR2_001_fastqc/Images/adapter_content.png' % self.sample)
                    plot.add_caption('R2 BEFORE trimming')
                with doc.create(SubFigure()) as plot:
                    if r2_summary_trim_dict.get('Adapter Content') == 'PASS':
                        colour = 'green'
                    elif r2_summary_trim_dict.get('Adapter Content') == 'WARN':
                        colour = 'orange'
                    elif r2_summary_trim_dict.get('Adapter Content') == 'FAIL':
                        colour = 'red'
                    else:
                        colour = 'black'
                    plot.add_image('%sR2_001.qfilter_fastqc/Images/adapter_content.png' % self.sample)
                    plot.add_caption(NoEscape(r'R2 AFTER trimming \textcolor{%s}{%s}'
                                              % (colour, r2_summary_trim_dict.get('Adapter Content'))))

        doc.append(Command('needspace', '1635em'))
        with doc.create(Section('BamStats')):
            with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))):
                doc.append(Command('centering'))
                with doc.create(SubFigure()) as plot:
                    plot.add_image('%s.bwa.drm.sorted.bam.stats-quals-hm.png' % self.sample)
                    plot.add_caption('Base quality per cycle')
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('hspace', '10 mm'))
                    plot.add_image('%s.bwa.drm.sorted.bam.stats-insert-size.png' % self.sample)
                    plot.add_caption('Fragment size')
                with doc.create(SubFigure()) as plot:
                    doc.append(Command('vspace', '10 mm'))
                    plot.add_image('%s.bwa.drm.sorted.bam.stats-quals2.png' % self.sample)
                    plot.add_caption('Quality per cycle')

        pdflatex = '/usr/local/texlive/2015/bin/x86_64-linux/pdflatex'
        doc.generate_pdf('%s' % self.sample, clean_tex=False, compiler=pdflatex)
        os.system('mv /home/cuser/PycharmProjects/amlpipeline/%s.pdf /media/sf_sarah_share/MiSeq_quality_outputs/'
                  % self.sample)

    def get_trimmed_data(self):
        # Get R1 trimmed
        r1_summary_trim_df = pd.read_table('%sR1_001.qfilter_fastqc/summary.txt' % self.sample, header=None,
                                           names=['Score', 'Parameter'], usecols=[0, 1])
        r1_scores_trim = r1_summary_trim_df['Score'].tolist()  # not currently used, may be needed
        r1_parameters_trim = r1_summary_trim_df['Parameter'].tolist()
        r1_summary_trim_dict = dict(zip(r1_parameters_trim, r1_scores_trim))
        r1_stats_trim_df = pd.read_table('%sR1_001.qfilter_fastqc/fastqc_data.txt' % self.sample, header=None,
                                         names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
        r1_properties_trim = r1_stats_trim_df['Property'].tolist()
        r1_values_trim = r1_stats_trim_df['Value'].tolist()
        r1_stats_trim_dict = dict(zip(r1_properties_trim, r1_values_trim))

        # Get R2 trimmed
        r2_summary_trim_df = pd.read_table('%sR2_001.qfilter_fastqc/summary.txt' % self.sample, header=None,
                                           names=['Score', 'Parameter'], usecols=[0, 1])
        r2_scores_trim = r2_summary_trim_df['Score'].tolist()  # not currently used, may be needed
        r2_parameters_trim = r2_summary_trim_df['Parameter'].tolist()
        r2_summary_trim_dict = dict(zip(r2_parameters_trim, r2_scores_trim))
        r2_stats_trim_df = pd.read_table('%sR2_001.qfilter_fastqc/fastqc_data.txt' % self.sample, header=None,
                                         names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
        r2_properties_trim = r2_stats_trim_df['Property'].tolist()
        r2_values_trim = r2_stats_trim_df['Value'].tolist()
        r2_stats_trim_dict = dict(zip(r2_properties_trim, r2_values_trim))

        return r1_summary_trim_dict, r1_stats_trim_dict, r2_summary_trim_dict, r2_stats_trim_dict

    def get_original_data(self):
        # Get R1
        r1_summary_df = pd.read_table('%sR1_001_fastqc/summary.txt' % self.sample, header=None,
                                      names=['Score', 'Parameter'], usecols=[0, 1])
        r1_scores = r1_summary_df['Score'].tolist()  # not currently used, may be needed
        r1_parameters = r1_summary_df['Parameter'].tolist()
        r1_summary_dict = dict(zip(r1_parameters, r1_scores))
        r1_stats_df = pd.read_table('%sR1_001_fastqc/fastqc_data.txt' % self.sample, header=None,
                                    names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
        r1_properties = r1_stats_df['Property'].tolist()
        r1_values = r1_stats_df['Value'].tolist()
        r1_stats_dict = dict(zip(r1_properties, r1_values))

        # Get R2
        r2_summary_df = pd.read_table('%sR2_001_fastqc/summary.txt' % self.sample, header=None,
                                      names=['Score', 'Parameter'], usecols=[0, 1])
        r2_scores = r2_summary_df['Score'].tolist()  # not currently used, may be needed
        r2_parameters = r2_summary_df['Parameter'].tolist()
        r2_summary_dict = dict(zip(r2_parameters, r2_scores))
        r2_stats_df = pd.read_table('%sR2_001_fastqc/fastqc_data.txt' % self.sample, header=None,
                                    names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
        r2_properties = r2_stats_df['Property'].tolist()
        r2_values = r2_stats_df['Value'].tolist()
        r2_stats_dict = dict(zip(r2_properties, r2_values))

        return r1_summary_dict, r1_stats_dict, r2_summary_dict, r2_stats_dict

pdf = CreateFastQCPDF('04-D15-22373-HT-Nextera-Myeloid-Val1-Repeat_S4_L001_')
pdf.create_pdf()
