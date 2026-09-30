import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, AllChem
import py3Dmol
import requests

st.set_page_config(page_title="Molecule Visualizer", page_icon="🧪", layout="wide")

st.markdown("""
    <style>
    .main { background-color: #0e1117; }
    .rule-pass { color: #10b981; }
    .rule-fail { color: #ef4444; }
    </style>
""", unsafe_allow_html=True)

st.title("🧪 Molecule Visualizer")
st.markdown("Enter a molecule name or a SMILES string to explore its structure and drug-likeness.")

EXAMPLES = {
    "Aspirin":    "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine":   "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen":  "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Dopamine":   "NCCc1ccc(O)c(O)c1",
    "Paracetamol":"CC(=O)Nc1ccc(O)cc1",
}

# ── Input area ──────────────────────────────────────────────────────────────
col_in, col_ex = st.columns([3, 1])

with col_in:
    mode = st.radio("Input type", ["Molecule name", "SMILES string"], horizontal=True)
    user_input = st.text_input(
        "Molecule" if mode == "SMILES string" else "Molecule name",
        placeholder="CC(=O)Oc1ccccc1C(=O)O" if mode == "SMILES string" else "aspirin",
    )

with col_ex:
    st.markdown("**Quick examples**")
    for name in EXAMPLES:
        if st.button(name, use_container_width=True):
            user_input = EXAMPLES[name]
            mode = "SMILES string"

# ── Resolve SMILES ───────────────────────────────────────────────────────────
def name_to_smiles(name: str):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES,CanonicalSMILES/JSON"
    try:
        r = requests.get(url, timeout=8)
    except Exception as e:
        st.error(f"Network error: {e}")
        return None
    if r.status_code != 200:
        st.error(f"Molecule '{name}' not found in PubChem. Check the spelling or use SMILES input.")
        return None
    props = r.json().get("PropertyTable", {}).get("Properties", [{}])[0]
    return props.get("IsomericSMILES") or props.get("CanonicalSMILES") or props.get("SMILES")

# ── 3D viewer ────────────────────────────────────────────────────────────────
def show_3d(mol, style="stick"):
    mol_h = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol_h, randomSeed=42)
    if result != 0:
        st.warning("Could not generate 3D coordinates for this molecule.")
        return
    AllChem.MMFFOptimizeMolecule(mol_h)
    sdf = Chem.MolToMolBlock(mol_h)

    view = py3Dmol.view(width=500, height=380)
    view.addModel(sdf, "mol")

    if style == "stick":
        view.setStyle({"stick": {}})
    elif style == "sphere":
        view.setStyle({"sphere": {"scale": 0.4}})
    elif style == "ball-stick":
        view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    elif style == "surface":
        view.setStyle({"stick": {}})
        view.addSurface(py3Dmol.VDW, {"opacity": 0.6, "colorscheme": "whiteCarbon"})

    view.setBackgroundColor("#0e1117")
    view.zoomTo()
    view.spin(True)

    html = view._make_html()
    components.html(html, height=400)

# ── Main logic ───────────────────────────────────────────────────────────────
if user_input:
    smiles = user_input.strip() if mode == "SMILES string" else name_to_smiles(user_input.strip())

    if not smiles:
        st.stop()

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        st.error("Invalid SMILES string — could not parse the molecule.")
        st.stop()

    st.divider()

    # ── Tabs: 2D / 3D / Properties ──────────────────────────────────────────
    tab2d, tab3d, tabprops = st.tabs(["🖼️ 2D Structure", "🔬 3D Viewer", "📊 Properties"])

    with tab2d:
        img = Draw.MolToImage(mol, size=(420, 320))
        col_img, col_smi = st.columns([2, 1])
        with col_img:
            st.image(img, use_container_width=True)
        with col_smi:
            st.markdown("**SMILES**")
            st.code(smiles, language=None)

    with tab3d:
        style = st.selectbox(
            "Display style",
            ["stick", "ball-stick", "sphere", "surface"],
            format_func=lambda s: {"stick": "Stick", "ball-stick": "Ball & Stick",
                                   "sphere": "Space-filling", "surface": "Surface"}[s],
        )
        show_3d(mol, style)
        st.caption("Drag to rotate · Scroll to zoom · The molecule spins automatically")

    with tabprops:
        mw    = Descriptors.MolWt(mol)
        logp  = Descriptors.MolLogP(mol)
        hbd   = rdMolDescriptors.CalcNumHBD(mol)
        hba   = rdMolDescriptors.CalcNumHBA(mol)
        tpsa  = rdMolDescriptors.CalcTPSA(mol)
        rb    = rdMolDescriptors.CalcNumRotatableBonds(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        heavy = mol.GetNumHeavyAtoms()

        c1, c2 = st.columns(2)
        c1.metric("Molecular Weight", f"{mw:.2f} g/mol")
        c2.metric("LogP (lipophilicity)", f"{logp:.2f}")
        c1.metric("H-Bond Donors", hbd)
        c2.metric("H-Bond Acceptors", hba)
        c1.metric("Rotatable Bonds", rb)
        c2.metric("TPSA", f"{tpsa:.1f} Å²")
        c1.metric("Rings", rings)
        c2.metric("Heavy Atoms", heavy)

        st.subheader("Lipinski Rule of Five")
        st.caption("Predicts whether a molecule could be orally bioavailable as a drug.")

        rules = [
            ("Molecular Weight ≤ 500 Da", mw <= 500),
            ("LogP ≤ 5",                  logp <= 5),
            ("H-Bond Donors ≤ 5",         hbd <= 5),
            ("H-Bond Acceptors ≤ 10",     hba <= 10),
        ]
        violations = sum(1 for _, ok in rules if not ok)

        for rule, ok in rules:
            st.markdown(f"{'✅' if ok else '❌'} {rule}")

        st.divider()
        if violations == 0:
            st.success("Drug-like — passes all Lipinski rules.")
        elif violations == 1:
            st.warning("1 violation — borderline drug-like.")
        else:
            st.error(f"{violations} violations — poor oral bioavailability predicted.")
