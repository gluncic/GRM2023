# GRMhor: A Tool for Automatic Annotation of Genome Monomer Higher Order Structure

**GRMhor** is a Python-based tool developed as part of the scientific paper titled "Efficient genome monomer higher order structure annotation and identification using the GRMhor algorithm". It provides functionalities for analyzing and visualizing the higher-order repeat (HOR) structure of genome monomers. This README provides an overview of the tool's functionality, usage instructions, and installation guidelines.

## Features

- **GRM Analysis**: Generates Genome Repeat Map (GRM) data from a sequence of monomers.
- **HOR Structure Visualization**: Visualizes the higher-order repeat (HOR) structure of monomers using a graphical representation.
- **Family Detection**: Groups monomers into families based on sequence similarity.
- **Family Merging**: Merges overlapping families of monomers to form larger, unified families.
- **File Input/Output**: Reads monomer data from a file and writes analysis results to output files.
- **Command-Line Interface**: Provides a command-line interface for easy access to functionality.

## Usage

1. **Installation**: Clone the repository and install the required dependencies using `pip`:

```bash
git clone https://github.com/yourusername/GRMhor.git
cd GRMhor
pip install -r requirements.txt
```

## Usage

1. **Running the Tool**: Execute the main script with the desired input file:

    ```bash
    python main_GRMhor.py input_file.txt --start 0 --pmax 60 --horpos
    ```

    Replace `input_file.txt` with the path to your input file. Adjust optional parameters `--start`, `--pmax`, and `--horpos` as needed:  
    **--start (default: 0):**  
      Defines the starting monomer in the sequence of monomers provided in input_file.txt. This is particularly useful for isolating smaller subsequences from long monomeric sequences,       simplifying the analysis.  
    **--pmax (default: 60):**  
      Specifies the maximum period displayed in the output diagrams, helping to clarify the visualization of HOR structures.  
    **--horpos (default: False):**  
      Prints the position (in base pairs) of the first monomer in each HOR unit, adding genomic context to the HOR structure.  


3. **Viewing Results**: Explore the input files and generated output files, including GRM diagrams, MD diagrams, and HOR structure visualizations, in the following directory: github.com/gluncic/GRM2023/tree/master/data.

## Citation

If you use **GRMhor** in your research work, please cite the corresponding paper: "Efficient genome monomer higher order structure annotation and identification using the GRMhor algorithm"

## License

This project is licensed under the [MIT License](LICENSE.md).

## Contribution

Contributions to **GRMhor** are welcome! Feel free to fork the repository, make improvements, and submit pull requests. Please ensure adherence to coding standards and maintain clear commit messages for better collaboration.

## Contact

For any questions, feedback, or support, please contact [matko.phy@pmf.hr](mailto:matko.phy@pmf.hr).
