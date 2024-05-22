# GRMhor: A Tool for Automatic Annotation of Genome Monomer Higher Order Structure

**GRMhor** is a Python-based tool developed as part of the scientific paper titled "GRMhor: a tool for automatic annotation genome monomer higher order structure". It provides functionalities for analyzing and visualizing the higher-order repeat (HOR) structure of genome monomers. This README provides an overview of the tool's functionality, usage instructions, and installation guidelines.

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
    python GRMhor.py input_file.txt --start 0 --pmax 60 --horpos
    ```

    Replace `input_file.txt` with the path to your input file. Adjust optional parameters `--start`, `--pmax`, and `--horpos` as needed.

2. **Viewing Results**: Check the generated output files, including GRM and HOR structure visualizations, in the specified directory.

## Citation

If you use **GRMhor** in your research work, please cite the corresponding paper:

## License

This project is licensed under the [MIT License](LICENSE.md).

## Contribution

Contributions to **GRMhor** are welcome! Feel free to fork the repository, make improvements, and submit pull requests. Please ensure adherence to coding standards and maintain clear commit messages for better collaboration.

## Contact

For any questions, feedback, or support, please contact [matko.phy@pmf.hr](mailto:matko.phy@pmf.hr).
