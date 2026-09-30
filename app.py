import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
import requests

st.set_page_config(page_title="Molecule Visualizer", page_icon="🧪", layout="wide")

st.markdown("""
    <style>
    .main { background-color: #0e1117; }
    .metric-box { background: #1e2130; border-radius: 10px; padding: 12px 16px; margin-bottom: 8px; }
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
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
    try:
        r = requests.get(url, timeout=8)
        if r.status_code == 200:
            return r.json()["PropertyTable"]["Properties"][0]["IsomericSMILES"]
        else:
            st.error(f"PubChem returned status {r.status_code} for '{name}'. Try SMILES input instead.")
            return None
    except Exception as e:
        st.error(f"Network error: {e}")
        return None

# ── Main logic ───────────────────────────────────────────────────────────────
if user_input:
    smiles = user_input.strip() if mode == "SMILES string" else name_to_smiles(user_input.strip())

    if not smiles:
        st.error("Molecule not found. Try a different name or switch to SMILES input.")
        st.stop()

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        st.error("Invalid SMILES string — could not parse the molecule.")
        st.stop()

    st.divider()
    left, right = st.columns([1, 1])

    # ── Structure ─────────────────────────────────────────────────────────
    with left:
        st.subheader("2D Structure")
        img = Draw.MolToImage(mol, size=(420, 320))
        st.image(img, use_container_width=True)
        st.code(smiles, language=None)

    # ── Properties ───────────────────────────────────────────────────────
    with right:
        st.subheader("Molecular Properties")

        mw   = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd  = rdMolDescriptors.CalcNumHBD(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        rb   = rdMolDescriptors.CalcNumRotatableBonds(mol)
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
            icon = "✅" if ok else "❌"
            st.markdown(f"{icon} {rule}")

        st.divider()
        if violations == 0:
            st.success("Drug-like — passes all Lipinski rules.")
        elif violations == 1:
            st.warning("1 violation — borderline drug-like.")
        else:
            st.error(f"{violations} violations — poor oral bioavailability predicted.")
