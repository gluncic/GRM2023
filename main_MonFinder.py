import sys
import time
from MonFinder import read_fasta_file, analyse_chromosome, find_monomers, join_direct_and_reverse_complement, write_monomers_file


def main():
    """
    The main function to execute the monomer finding process.
    """
    start_time = time.time()

    if len(sys.argv) != 2:
        print("Usage: python3 MonFinder.py path_to_file")
        sys.exit(1)

    file_name = sys.argv[1]
    monomer_cons = "TCAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAGTTGAACCTTCCTTTTGATAGAGCAGTTTTGAAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGAGCGCTTTGAGGCCTTCGGTGGAAAAGGAAATATCTTCACATAAAAACTAGACAGAAGCATTC"

    name, seq = read_fasta_file(file_name)
    analyse_chromosome(seq)

    print("Searching for monomers in direct orientation")
    monomers_d = find_monomers(seq, monomer_cons, 'd')
    print("Searching for monomers in reverse complement orientation")
    monomers_r = find_monomers(seq, monomer_cons, 'r')
    monomers = join_direct_and_reverse_complement(monomers_d, monomers_r)
    write_monomers_file(f"{file_name.split('.', 1)[0]}.mon", monomers)

    print(f"--- {time.time() - start_time} seconds ---")


if __name__ == "__main__":
    main()
