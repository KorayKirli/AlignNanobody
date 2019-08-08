from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
import glob
import sys
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import numpy as np
import panel as pn
pn.extension()


def get_files(folder_name, ext):
    full_folder = '../data/' + folder_name
    file_list = glob.glob(full_folder + '/*.' + ext)
    num_files = len(file_list)
    if num_files:
        print(len(file_list), 'found in data folder', folder_name)
        return file_list
    else:
        raise Exception('no files were found in {} with extension {}'.format(
            folder_name, ext))


def extract_protein(dna_file, min_pro_len, only_starting_with_M, print_each):
    records = SeqIO.parse(dna_file, 'fasta')
    f_name = dna_file.split('/')[-1].split('.')[0]
    proteins = []
    for record in records:
        for strand, seq in (1, record.seq), (-1,
                                             record.seq.reverse_complement()):
            for frame in range(3):
                length = 3 * ((len(seq) - frame) // 3)
                for pro in seq[frame:frame + length].translate(
                        table=1).split("*")[:-1]:
                    if only_starting_with_M:
                        pro = pro[pro.find('M'):]
                    if len(pro) > min_pro_len:
                        proteins.append(pro)

#                     REPORTING
#                     pos = seq[frame:frame+length].translate(table=1).find(orf)*3 + frame +1
#                         print("{}...{} - length {}, strand {}, frame {}, pos {}, name {}".format\
#                            (orf, orf[:3], orf[-3:], len(orf)*3+3, strand, frame, pos, record.id))
    if not proteins:
        if print_each:
            print(f_name, 'protein size: 0', ', No protein found')
        return
    else:
        longest = max(proteins, key=len)
        xno = longest.count('X')
        if print_each:
            print(f_name, 'protein size:', len(longest), ', number of X:', xno)
        return {'fname': f_name, 'seq': longest, 'Xno': xno}


def store_proteins(project_name, cases):
    f_name = '../data/' + project_name + '.fasta'
    with open(f_name, 'w') as out_file:
        for a_case in cases:
            out_file.write('>' + a_case['fname'] + '\n')
            out_file.write(str(a_case['seq']) + '\n')
    return (f_name)


def run_msa(protein_file):
    muscle_exe = ""
    platform = sys.platform
    exe_map = {
        'darwin': './muscle/muscle3.8.31_i86darwin64',
        'linux': './muscle/muscle3.8.31_i86linux64',
        'win32': './muscle/muscle3.8.31_i86win32.exe',
        'cygwin': './muscle/muscle3.8.31_i86cygwin32.exe'
    }
    muscle_exe = exe_map.get(platform)
    if not muscle_exe:
        raise Exception('no muscle exe can be located')

    out_file = protein_file.replace('.fasta', '_align.fasta')
    muscle_cline = MuscleCommandline(
        muscle_exe, input=protein_file, out=out_file)
    stdout, stderr = muscle_cline()
    return out_file


def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    clrs = {}
    mapper = {
        'GAST': 'orange',
        'DE': 'red',
        'KR': 'blue',
        'NQR': 'magenta',
        'CVILPFYMW': 'green',
        'X': 'gray'}
    clrs = {i: mapper[comb] for comb in mapper for i in comb}
    colors = [clrs.get(i, 'white') for i in text]
    return colors


def view_alignment(aln, fontsize="9pt", plot_width=1500):
    """Bokeh sequence alignment view
    Modified from http://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner"""
    # make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)
    N = len(seqs[0])
    S = len(seqs)
    width = .4

    x = np.arange(1, N + 1)
    y = np.arange(0, S, 1)
    # creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    # flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    # use recty for rect coords with an offset
    recty = gy + .5
    h = 1 / S
    # now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(
        dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs) * 15 + 50
    x_range = Range1d(0, N + 1, bounds='auto')
    viewlen = N + 1
    if N > 250:
        viewlen = 250
    else:
        viewlen = N
    # view_range is for the close up view
    view_range = (0, viewlen)
    tools = "xpan, xwheel_zoom, reset, save"

    # entire sequence view (no text, with zoom)
    p = figure(
        title=None,
        plot_width=plot_width,
        plot_height=50,
        x_range=x_range,
        y_range=(0, S),
        tools=tools,
        min_border=0,
        toolbar_location='below')
    rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color=None,
        fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False

    # sequence text view with ability to scroll along x axis
    p1 = figure(
        title=None,
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=view_range,
        y_range=ids,
        tools="xpan,reset",
        min_border=0,
        toolbar_location='below')   # lod_factor=1)
    glyph = Text(
        x="x",
        y="y",
        text="text",
        text_align='center',
        text_color="black",
        text_font="monospace",
        text_font_size=fontsize)
    rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color=None,
        fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = gridplot([[p], [p1]], toolbar_location='below')
    return p
