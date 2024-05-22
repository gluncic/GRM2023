import edlib
from tkinter import Canvas, Tk
import seaborn as sns
import colorcet as cc
import matplotlib.pyplot as plt




class Monomer:
    """
    Represents a monomer with properties like position, distance to next, divergence scores, orientation, and sequence.
    Additional properties include family, row, column positions, and monomer number for internal tracking.
    """

    def __init__(self, pos, dst, div, div2, ort, seq):
        self.pos = int(pos)
        self.dst = int(dst)
        self.div = float(div)
        self.div2 = int(div2)
        self.ort = ort
        self.seq = seq
        self.family = []
        self.row = -1
        self.col = -1
        self.colStart = -1
        self.monNo = -1

    def __str__(self):
        return f"{self.pos} {self.dst} {self.div:.2f} {self.seq}"


def draw_grm_and_mdd(series, monomers, s_file_name, b_block_lines, xmax, ymax, xtics_period, ytics_period):
    """
    Draws Genome Repeat Map (GRM) and Monomer Distance Distribution (MDD) based on the analysis of monomer sequences.

    Parameters:
    - series: List[int] - The sequence of monomer indices for analysis.
    - monomers: List[Monomer] - The list of Monomer objects.
    - s_file_name: str - The base name for output files.
    - b_block_lines: bool - Flag to draw block lines.
    - xmax: int - Maximum x-axis value for plots.
    - ymax: int - Maximum y-axis value for plots.
    - xtics_period: int - X-axis tick period.
    - ytics_period: int - Y-axis tick period.
    """
    freq, frag, gap, gap_len = grm(series, monomers)
    sorted_freq = sorted(freq, reverse=True)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))
    plt.rcParams['font.serif'] = 'Helvetica'
    plt.subplots_adjust(hspace=0.23)
    ax1.tick_params(labelsize=30)
    ax2.tick_params(labelsize=30)

    ax1.plot(freq, color='black')
    for i in range(7):
        idx = freq.index(sorted_freq[i])
        ax1.text(idx, sorted_freq[i], idx, ha='center', va='bottom', fontsize=20)

    ax1.set_xlabel('Period', fontsize=30)
    ax1.set_ylabel('Frequency', fontsize=30)
    ax1.set_xlim(0, xmax)
    ax1.set_ylim(0, max(freq) + 10 * max(freq) / 100)
    ax1.set_xticks(range(0, xmax + 1, 5))

    ax2.plot(frag, marker='o', fillstyle='full', markersize=0.4, linestyle='', color='black')
    if b_block_lines:
        for i in range(len(gap)):
            ax2.plot([gap[i], gap[i]], [0, ymax], color='red', linewidth=0.1)
            ax2.text(gap[i], ymax, gap_len[i], ha='center', va='bottom', fontsize=2)

    ax2.set_xlabel('Index', fontsize=30)
    ax2.set_ylabel('Period', fontsize=30)
    ax2.grid(True, which='major', linestyle='dotted', linewidth=0.3, color='gray')

    ax2.set_xlim(0, len(frag))
    ax2.set_ylim(0, ymax)
    ax2.set_yticks(range(0, ymax + 1, 10))

    plt.savefig(f"{s_file_name}.GRM_MDD.pdf", format='pdf', bbox_inches='tight', pad_inches=0.01)


def grm(series, monomers):
    """
    Generates Genome Repeat Map (GRM) data from a series of monomer sequences.

    Parameters:
    - series: List[int] - The sequence of monomer indices for analysis.
    - monomers: List[Monomer] - The list of Monomer objects.

    Returns:
    - Tuple: A tuple containing frequency, fragment distances, gap positions, and gap lengths.
    """
    step = 1
    freq = [0] * 61
    frag = [0] * len(series)
    for i in range(len(series) - step - 1):
        match = series[i:i + step]
        for j in range(i + step, len(series) - step):
            subject = series[j:j + step]
            if match == subject:
                fragment = j - i
                frag[i] = fragment
                if fragment < 61:
                    freq[fragment] += 1
                break
    gap = []
    gap_len = []
    for i in range(len(monomers) - 1):
        if monomers[i + 1].pos - monomers[i].pos > 1000:
            gap.append(i + 1)
            gap_len.append(monomers[i + 1].pos - monomers[i].pos)

    return freq, frag, gap, gap_len


def draw_hor_structure(monomers, file_name, b_numbers, b_position_marks_blocks, b_mers_marks, b_alpha_positions,
                       f_cube_proportions):
    """
    Visualizes the Higher Order Repeat (HOR) structure of monomers using a graphical representation.

    Parameters:
    - monomers: List[Monomer] - The list of Monomer objects.
    - file_name: str - The base name for the output file.
    - b_numbers: bool - Flag to include monomer sequence numbers.
    - b_position_marks_blocks: bool - Flag to mark positions with blocks.
    - b_mers_marks: bool - Flag to include tags at the first occurrence above each monomer.
    - b_alpha_positions: bool - Flag to include position of the first monomer in the HOR.
    - f_cube_proportions: float - Factor to adjust the size of each square in the visualization.

    Returns:
    - List[int] - The series of monomer indices used for the GRM analysis.
    """
    x, y = [], []
    z = {}
    bx, by, mx, my = [], [], [], []
    pos_y, pos_pos = [], []
    j = i = 10
    pos_y.append(i)
    pos_pos.append(monomers[0].pos)
    last = maxlast = -1
    fill_and_fit_columns(monomers)
    for alpha in monomers:
        x.append(alpha.col)

        if alpha.col <= last:
            i += 1
            pos_y.append(i)
            pos_pos.append(alpha.pos)

        if alpha.col > maxlast and len(alpha.family) > 1:
            maxlast = alpha.col
            mx.append(alpha.col)
            my.append(i)
        last = alpha.col
        y.append(i)
        z[(alpha.col, i)] = j
        j += 1

        if alpha.dst > 1000:
            bx.append(alpha.col)
            by.append(i)

    palette = sns.color_palette(cc.palette.glasbey_light, n_colors=10000).as_hex()
    c = [i if x.count(i) > 1 else 0 for i in x]

    w = h = 0.8 * f_cube_proportions
    shift = 10 * f_cube_proportions
    app = Tk()
    canvas = Canvas(app, width=max(x) * f_cube_proportions + 10 + shift, height=max(y) * f_cube_proportions + 20)
    canvas.pack()

    for i in range(len(x)):
        canvas.create_rectangle(x[i] * f_cube_proportions + shift, y[i] * f_cube_proportions,
                                x[i] * f_cube_proportions + w + shift, y[i] * f_cube_proportions + h, outline='black',
                                fill=palette[c[i]], width=0.07 * f_cube_proportions)
        if b_numbers:
            canvas.create_text(x[i] * f_cube_proportions + shift, y[i] * f_cube_proportions, text=i,
                               font=("Helvetica", 1))

    for i in range(len(bx)):
        canvas.create_oval(bx[i] * f_cube_proportions + shift, by[i] * f_cube_proportions,
                           bx[i] * f_cube_proportions + w + shift, by[i] * f_cube_proportions + h, fill='white',
                           outline='red', width=0.07 * f_cube_proportions)
        if b_position_marks_blocks:
            txt = f"{z[bx[i], by[i]]} ({monomers[z[bx[i], by[i]]].pos})"
            canvas.create_text(bx[i] * f_cube_proportions + shift, by[i] * f_cube_proportions, text=txt,
                               font=("Helvetica", 1))
            if z[bx[i], by[i]] > 0:
                txt = f"{z[bx[i], by[i]] - 1} ({monomers[z[bx[i], by[i]] - 1].pos})"
                canvas.create_text(bx[i] * f_cube_proportions + shift, (by[i] - 1) * f_cube_proportions, text=txt,
                                   font=("Helvetica", 1))

    if b_alpha_positions:
        for i, pos in enumerate(pos_y):
            canvas.create_text(shift - 1 * f_cube_proportions, pos * f_cube_proportions + h / 2, text=f"{pos_pos[i]:,}",
                               font=("Helvetica", 1))

    if b_mers_marks:
        for i, mx_val in enumerate(mx):
            canvas.create_text(mx_val * f_cube_proportions + shift + w / 2, my[i] * f_cube_proportions - 1,
                               text=f"m{i + 1}", font=("Helvetica", 1))

    canvas.update()
    canvas.postscript(file=f"{file_name}.HORscheme.ps", colormode='color',
                      width=max(x) * f_cube_proportions + 10 + shift, height=max(y) * f_cube_proportions + 20)

    return x


def fill_and_fit_columns(monomers):
    """
    Assigns columns to monomers based on their families for visualization purposes.

    Parameters:
    - monomers: List[Monomer] - The list of Monomer objects.
    """
    i = 0
    for mon in monomers:
        if mon.family:
            for mem in mon.family:
                monomers[mem].col = monomers[mem].colStart = i
            i += 1
    if monomers[-1].col < monomers[0].col:
        monomers[-1].col = i


def find_families(monomers, limit):
    """
    Groups monomers into families based on sequence similarity.

    Parameters:
    - monomers: List[Monomer] - The list of Monomer objects to analyze.
    - limit: int - The maximum divergence percentage allowed for monomers to be considered part of the same family.
    """
    no = len(monomers)
    progress_step = max(no // 20, 1)
    for i in range(no - 1):
        if i % progress_step == 0:
            progress_percentage = i / no * 100
            print(f"Progress: {progress_percentage:5.1f}%")
        for j in range(i, no):
            ed = edlib.align(monomers[i].seq, monomers[j].seq)
            edp = ed["editDistance"] / len(monomers[i].seq) * 100
            if edp < limit:
                monomers[i].family.append(j)


def join_families_v03(monomers):
    """
    Merges overlapping families of monomers to form larger, unified families.

    Parameters:
    - monomers: List[Monomer] - The list of Monomer objects to analyze and merge families.
    """
    for i in range(len(monomers)):
        if len(monomers[i].family) > 1:
            for j in range(len(monomers)):
                if i == j:
                    continue
                if set(monomers[i].family) & set(monomers[j].family):
                    union = sorted(set(monomers[i].family) | set(monomers[j].family))
                    monomers[i].family = monomers[j].family = []
                    monomers[union[0]].family = union


def read_monomers_file(file_name, start):
    """
    Reads a file containing monomer data and returns a list of Monomer objects.

    Parameters:
    - file_name: str - The path to the file to read.
    - start: int - The starting index to begin processing from the file.

    Returns:
    - List[Monomer] - A list of Monomer objects read from the file.
    """
    monomers = []
    with open(file_name) as file:
        for i, line in enumerate(file.readlines()):
            if i >= start:
                parts = line.split()
                if len(parts) == 6:
                    monomers.append(Monomer(*parts))
                else:
                    monomers.append(Monomer(i, len(parts[0]), 0, 0, 'd', parts[0].strip()))
    return monomers
