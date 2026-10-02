# 🔬 Molecule Visualizer

An interactive web app for visualizing molecular structures and analyzing drug-likeness properties.

Built with **RDKit** and **Streamlit** as part of my bioinformatics & cheminformatics learning journey.

## Features

- 🔍 Search by **molecule name** (looks up PubChem) or paste a **SMILES string** directly
- 🖼️ **2D structure** rendering using RDKit
- 🔬 **Interactive 3D viewer** (py3Dmol) — drag to rotate, scroll to zoom, auto-spin
  - 4 display styles: Stick, Ball & Stick, Space-filling, Surface
- 📊 **Molecular properties**: molecular weight, LogP, H-bond donors/acceptors, TPSA, rotatable bonds, ring count
- 🎯 **Drug-likeness tab**:
  - **QED score** (0–1 gauge) — Quantitative Estimate of Drug-likeness
  - **Property radar chart** — all 6 properties plotted against their Lipinski/Veber thresholds
  - **Rules scorecard** — Lipinski Rule of Five + Veber rules with pass/fail and actual values
- ⚡ Quick-load examples: Aspirin, Caffeine, Ibuprofen, Dopamine, Paracetamol

## Getting Started

```bash
git clone https://github.com/marishachem/molecule-visualizer.git
cd molecule-visualizer
pip install -r requirements.txt
python3 -m streamlit run app.py
```

Opens at `http://localhost:8501`

## Requirements

- Python 3.8+
- RDKit, Streamlit, Requests, Pillow, py3Dmol, matplotlib, numpy

## What I Learned

- How SMILES notation represents molecular structures
- Using RDKit for cheminformatics (2D drawing, 3D conformer generation, descriptors, QED)
- Generating 3D coordinates with `AllChem.EmbedMolecule` + MMFF force field optimization
- Interactive 3D visualization with py3Dmol embedded in Streamlit
- Lipinski's Rule of Five and Veber's rules in drug discovery
- Radar/spider charts with matplotlib for multi-property visualization
- Integrating the PubChem REST API

## Next Steps

- [ ] Download molecule as PNG/SDF
- [ ] Batch screening of multiple molecules
- [ ] Side-by-side molecule comparison
