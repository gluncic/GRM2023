
import edlib


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


def join_direct_and_reverse_complement(monomers1, monomers2):
    """
    Merges two lists of Monomer objects, sorting them by position and removing duplicates based on orientation rules.

    Parameters:
    - monomers1: List[Monomer] - The first list of monomers.
    - monomers2: List[Monomer] - The second list of monomers, typically the reverse complement monomers.

    Returns:
    - List[Monomer] - The merged and processed list of monomers.
    """
    monomers1.extend(monomers2)
    monomers1.sort(key=lambda x: x.pos)

    monomers = []
    skip = False
    for i in range(len(monomers1) - 2):
        if skip:
            skip = False
        else:
            monomers.append(monomers1[i])
            if monomers1[i].ort != monomers1[i + 1].ort and monomers1[i].ort == monomers1[i + 2].ort:
                skip = True

    monomers.append(monomers1[-2])
    monomers.append(monomers1[-1])
    set_distances(monomers)
    return monomers


def find_monomers(seq, monomer_cons, ort):
    """
    Identifies and returns a list of Monomer objects found in a given sequence.

    Parameters:
    - seq: str - The sequence in which to find monomers.
    - monomer_cons: str - The consensus sequence of the monomer.
    - ort: str - The orientation ('d' for direct, 'r' for reverse complement).

    Returns:
    - List[Monomer] - A list of identified Monomer objects.
    """
    monomers = []
    if ort == "r":
        seq = revcom(seq)

    progress_step = max((len(seq) - len(monomer_cons) - 1) // 20, 1)
    for i in range(len(seq) - len(monomer_cons) - 1):
        if i % progress_step == 0:
            progress_percentage = i / (len(seq) - len(monomer_cons) - 1) * 100
            print(f"Progress: {progress_percentage:5.1f}%")
        ed = edlib.align(monomer_cons, seq[i: i + len(monomer_cons)])
        edp = ed["editDistance"] / len(monomer_cons) * 100
        if edp < 30:
            ed = edlib.align(monomer_cons[:10], seq[i: i + 10])
            monomers.append(Monomer(i, 0, edp, ed["editDistance"], ort, seq[i: i + 10]))
    monomers = find_min_alphas(monomers)
    set_distances(monomers)
    monomers = remove_small_distances(monomers)
    set_distances(monomers)
    set_sequences(monomers, seq)
    if ort == "r":
        set_back_rc_positions(monomers, len(seq))
        monomers.sort(key=lambda x: x.pos)
    return monomers


def set_back_rc_positions(monomers, seq_len):
    """
    Adjusts the positions of monomers to their original locations in the reverse complement sequence.

    Parameters:
    - monomers: List[Monomer] - The list of monomers to adjust.
    - seq_len: int - The length of the original sequence.
    """
    for monomer in monomers:
        monomer.pos = seq_len - monomer.pos


def find_min_alphas(monomers):
    """
    Filters a list of monomers to include only those with minimum divergence scores compared to their immediate neighbors.

    Parameters:
    - monomers: List[Monomer] - The original list of monomers.

    Returns:
    - List[Monomer] - The filtered list of monomers.
    """
    monomers_min = []
    for i in range(1, len(monomers) - 1):
        if monomers[i].div <= monomers[i - 1].div and monomers[i].div <= monomers[i + 1].div:
            monomers_min.append(monomers[i])
    return monomers_min


def set_sequences(monomers, seq):
    """
    Sets the sequence for each monomer based on its position and the distance to the next monomer.

    Parameters:
    - monomers: List[Monomer] - The list of monomers to update.
    - seq: str - The original sequence from which the monomers were identified.
    """
    for i in range(len(monomers) - 1):
        if monomers[i + 1].dst < 180:
            monomers[i].seq = seq[monomers[i].pos: monomers[i + 1].pos]
        else:
            monomers[i].seq = seq[monomers[i].pos: monomers[i].pos + 171]
    if monomers:
        monomers[-1].seq = seq[monomers[-1].pos: monomers[-1].pos + 171]


def remove_small_distances(monomers):
    """
    Filters out monomers that are too close to each other, based on a distance threshold.

    Parameters:
    - monomers: List[Monomer] - The list of monomers to filter.

    Returns:
    - List[Monomer] - The filtered list of monomers.
    """
    groups = []
    group = []
    for monomer in monomers:
        if monomer.dst > 20:
            if group: groups.append(group)
            group = []
        group.append(monomer)
    if group: groups.append(group)
    mons = [find_alpha_smallest_div10(g) for g in groups if g]
    return mons


def find_alpha_smallest_div10(group):
    """
    Finds the monomer with the smallest divergence score in a group of monomers.

    Parameters:
    - group: List[Monomer] - The group of monomers to evaluate.

    Returns:
    - Monomer - The monomer with the smallest divergence score.
    """
    best_score = 10
    best_monomer = group[0]  # Default to first monomer if no better found
    for monomer in group:
        if monomer.div2 <= best_score:
            best_monomer = monomer
            best_score = monomer.div2
    return best_monomer


def set_distances(monomers):
    """
    Calculates and sets the distance to the next monomer for each monomer in a list.

    Parameters:
    - monomers: List[Monomer] - The list of monomers to update.
    """
    last_pos = 0
    for monomer in monomers:
        monomer.dst = monomer.pos - last_pos
        last_pos = monomer.pos


def read_fasta_file(file_name):
    """
    Reads a FASTA file and returns the sequence title and the sequence itself.

    Parameters:
    - file_name: str - The path to the FASTA file.

    Returns:
    - tuple: A tuple containing the sequence title and the sequence.
    """
    with open(file_name, 'r') as file:
        title = file.readline().strip()
        name = title.split('|')[3] if '|' in title else title
        genome = file.read().replace('\r', '').replace('\n', '')
    return name, genome


def write_monomers_file(file_name, monomers):
    """
    Writes a list of monomers to a file.

    Parameters:
    - file_name: str - The path to the output file.
    - monomers: List[Monomer] - The list of monomers to write.
    """
    if monomers:
        with open(file_name, "w") as fp:
            for mon in monomers:
                fp.write(f"{mon.pos} {mon.dst} {mon.div:.2f} {mon.div2} {mon.ort} {mon.seq}\n")


def complement(base):
    """
    Returns the complement of a nucleotide base.

    Parameters:
    - base: str - The nucleotide base.

    Returns:
    - str - The complement of the base.
    """
    return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'X': 'X'}.get(base, 'N')


def revcom(s):
    """
    Returns the reverse complement of a DNA sequence.

    Parameters:
    - s: str - The DNA sequence.

    Returns:
    - str - The reverse complement of the sequence.
    """
    return "".join(complement(base) for base in reversed(s))


def analyse_chromosome(seq):
    """
    Prints the count of each nucleotide in a given sequence.

    Parameters:
    - seq: str - The DNA sequence to analyze.
    """
    print('A:', seq.count('A'))
    print('C:', seq.count('C'))
    print('G:', seq.count('G'))
    print('T:', seq.count('T'))
    print('N:', seq.count('N'))
    print("Sequence Length:", len(seq))

