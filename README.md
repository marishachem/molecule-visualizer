# 🧪 Molecule Visualizer

An interactive web app for visualizing molecular structures and analyzing drug-likeness properties.

Built with **RDKit** and **Streamlit** as part of my bioinformatics & cheminformatics learning journey.

## Features

- 🔍 Search by **molecule name** (looks up PubChem) or paste a **SMILES string** directly
- 🖼️ Renders **2D molecular structure** using RDKit
- 📊 Calculates key **molecular properties**:
  - Molecular weight
  - LogP (lipophilicity)
  - H-bond donors & acceptors
  - TPSA, rotatable bonds, ring count
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
- Using RDKit for cheminformatics (drawing, descriptors)
- Lipinski's Rule of Five and its role in drug discovery
- Integrating the PubChem REST API
- Building interactive data apps with Streamlit

## Next Steps

- [ ] Add 3D structure viewer (py3Dmol)
- [ ] Download molecule as PNG/SDF
- [ ] Batch screening of multiple molecules
- [ ] Drug-likeness score visualization
