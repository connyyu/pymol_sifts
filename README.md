# Automatic protein labelling in PyMOL.

This repository contains a PyMOL plugin that automatically label proteins chains by their UniProt IDs.

The tool retreives the SIFTS mapping (jointly provided by UniProt and PDBe) using the PDBe API. It automatically colours each protein chain and creates a named selection using the UniProt ID of the protein.

## Features

- Each protein is clearly coloured and labelled!
- Limitation: not compatible with chimeric and non-protein constructs.

## Installation
In PyMOL 3.0+ or Open-Source PyMOL, go to _Plugin_ -> _Plugin Manager_ <br>
In the _Install New Plugin_ tab, under _Install from PyMOLWiki or any URL:_

Fetch this URL:<br>
[https://github.com/connyyu/list_proteins/blob/main/sifts.py
](https://github.com/connyyu/pymol_sifts/blob/main/sifts.py)

### Usage
`sifts <object_name> [<object_name2> ...]`

<object_name> should match with PDB identifier.

<img src="https://github.com/user-attachments/assets/85e7c8d8-86d3-4378-9bdd-bf9bd85706d2" alt="SIFTS PyMOL plugin" width="600"/>

## Compatibility

Tested and fully functional on:
- PyMOL Version 3.1.6.1
- OpenSource PyMOL Version 3.2.0a

## Disclaimer & Attribution

PyMOL is a trademark of Schrödinger, LLC. This is an independent community project and is not officially affiliated with or endorsed by Schrödinger.

This tool leverages the SIFTS mapping (jointly provided by UniProt and PDBe)<br>
More about SIFTS: https://www.ebi.ac.uk/pdbe/docs/sifts/

More about the example protein structures:<br>
[Structure of human PINK1 at a mitochondrial TOM-VDAC array.](https://doi.org/10.1126/science.adu6445)
    
## Author

- **Conny Yu** – [GitHub Profile](https://github.com/connyyu)  
  _Mar 2026_
