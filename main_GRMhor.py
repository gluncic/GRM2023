import argparse
import time
from GRMhor import read_monomers_file,find_families,join_families_v03,draw_hor_structure,draw_grm_and_mdd
def main(input_file, start=None, pmax=None, horpos=None):
    """
    The main function to process the input file and generate visualization of monomers and their structures.

    Parameters:
    - input_file: str - Path to the input file containing monomer data.
    - start: int - Starting index for processing monomers.
    - pmax: int - Maximum period to be considered for analysis.
    - horpos: bool - Flag to include position of the first monomer in the HOR in the output.
    """
    start_time = time.time()

    monomers = read_monomers_file(input_file, start)
    print(f"{input_file} -> No monomers = {len(monomers)}")

    find_families(monomers, 5)
    join_families_v03(monomers)

    series = draw_hor_structure(monomers, input_file, b_numbers=False, b_position_marks_blocks=False,
                                b_mers_marks=False, b_alpha_positions=horpos, f_cube_proportions=1)
    draw_grm_and_mdd(series, monomers, input_file, b_block_lines=False, xmax=pmax, ymax=pmax, xtics_period=2000,
                     ytics_period=5)

    print(f"--- {time.time() - start_time} seconds ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='GRMhor: a tool for automatic annotation of genome monomer higher order structure')
    parser.add_argument('input_file', help='Required input file name')
    parser.add_argument('--start', type=int, default=0, help='Starting monomer in the sequence (default=0)')
    parser.add_argument('--pmax', type=int, default=60, help='Maximum value of the displayed period')
    parser.add_argument('--horpos', action='store_true', default=False,
                        help='Prints the position of the first monomer in the HOR')
    args = parser.parse_args()

    main(args.input_file, args.start, args.pmax, args.horpos)

