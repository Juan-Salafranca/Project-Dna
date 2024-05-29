# DNA and protein analysis tool : DNAPROT

<div style="text-align: center;">
  <img src="./assets/project_logo.jpg" alt="Alt text" width="400" height="400">
</div>

## Introduction
This project was made in the context of the CH-200 Practical Programming for chemistry course at EPFL by two second year Bachelor students in Chemistry and Chemical engineering. 

## Authors
- Marie Van Rossum, Bsc in Chemistry and chemical engeneering at EPFL. [<img src="./assets/GitHubSymb.png" alt="Profile Picture" width="60">](https://github.com/Marie-vanrossum)
- Juan Salafranca Martinez, Bsc in Chemistry and chemical engeneering at EPFL. [<img src="./assets/GitHubSymb.png" alt="Profile Picture" width="60">](https://github.com/Juan-Salafranca)

## Overview
The Protein Analysis Package is a Python package designed to perform various analyses on ADN and their respective protein sequences. The package includes functions to convert ADN sequences into protein sequences by evaluating all possible frames of ADN. It is also able to analyse the protein sequences and calculate the hydrophobicity score, molecular weight, likelihood of secondary structure configurations (beta-sheet, alpha-helix, and beta-turn), and the retention coefficient in High-Performance Liquid Chromatography (HPLC) based on the amino acid sequences.

## Features

- Shine-Dalgarno identification sequence: Identify Shine-Dalgarno sequence among a long sequence of ADN and give all possible protein convertible sections.
- ADN to Protein conversion: Convert the ADN into the protein sequences using different reading frames.
- Hydrophobicity Score Calculation: Compute the hydrophobicity score of a protein sequence using Kyte & Doolittle's scale.
- Molecular Weight Calculation: Determine the molecular weight of a protein sequence based on the molecular weights of individual amino acids.
- Secondary Structure Likelihood Calculation: Evaluate the likelihood of a protein sequence forming beta-sheets, alpha-helices, or beta-turn using the Chou and Fasman techniques
- HPLC Retention Coefficient Calculation: Calculate the retention coefficient of a protein sequence in HPLC using given retention values for amino acids.
- Polarity evaluation: Polarity score calculation based on the Zimmerman scale.
 
## Installation

```
pip install DNAPROT
```

The package uses the following dependencies:

-Pandas (2.2.2)

-Openpyxl (3.1.3)

## Usage

Attention: for the transcription functions to work, a file containing the genetic code file must be in the same repo as the script. One can use the file provided alongside this package or use a different code for other genetic codes (ex: mitochondrial ADN code, etc...)

Start by importing the package:

```
import DNAPROT_analysis as DNAPROT
```



## License
This project is open-source and released under the [MIT License](./LICENSE.txt).

## Sources
This project used different science papers to achieve its goal and determine properties of proteins. These sources can be found in the document: [Sources](./SOURCES.md).

## Contact
For any questions or issues, please contact juan.salafrancamartinez@epfl.ch or marie.vanrossum@epfl.ch
