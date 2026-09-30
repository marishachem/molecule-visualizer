# 🧪 Molecule Visualizer

An interactive web app for visualizing molecular structures and analyzing drug-likeness properties.

Built with **RDKit** and **Streamlit** as part of my bioinformatics & cheminformatics learning journey.

## Features

- 🔍 Search by **molecule name** (looks up PubChem) or paste a **SMILES string** directly
- 🖼️ **2D structure** rendering using RDKit
- 🔬 **Interactive 3D viewer** (py3Dmol) — drag to rotate, scroll to zoom, auto-spin
  - 4 display styles: Stick, Ball & Stick, Space-filling, Surface
- 📊 Calculates key **molecular properties**:
  - Molecular weight, LogP (lipophilicity)
  - H-bond donors & acceptors, TPSA, rotatable bonds, ring count
- 💊 **Lipinski Rule of Five** checker — predicts oral drug-likeness
- ⚡ Quick-load examples: Aspirin, Caffeine, Ibuprofen, Dopamine, Paracetamol

## Demo

| Input | Output |
|-------|--------|
| Type `aspirin` → | 2D structure + properties + Lipinski check |
| Paste SMILES `CC(C)=O` → | Same for acetone |

## Getting Started

```bash
git clone https://github.com/marishachem/molecule-visualizer.git
cd molecule-visualizer
pip install -r requirements.txt
python3 -m streamlit run app.py
```

## Requirements

- Python 3.8+
- RDKit
- Streamlit
- Requests
- Pillow

## What I Learned

- How SMILES notation represents molecular structures
- Using RDKit for cheminformatics (2D drawing, 3D conformer generation, descriptors)
- Generating 3D coordinates with `AllChem.EmbedMolecule` + MMFF force field optimization
- Interactive 3D visualization with py3Dmol embedded in Streamlit
- Lipinski's Rule of Five and its role in drug discovery
- Integrating the PubChem REST API

## Next Steps

- [ ] Download molecule as PNG/SDF
- [ ] Batch screening of multiple molecules
- [ ] Drug-likeness score visualization
