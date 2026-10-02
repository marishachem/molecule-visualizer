import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, AllChem, QED
import py3Dmol
import requests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

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

    # ── Tabs: 2D / 3D / Properties / Drug-likeness ──────────────────────────
    tab2d, tab3d, tabprops, tabdrug = st.tabs(["🖼️ 2D Structure", "🔬 3D Viewer", "📊 Properties", "🎯 Drug-likeness"])

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

    with tabdrug:
        qed = QED.qed(mol)

        # ── QED gauge ───────────────────────────────────────────────────────
        st.subheader("QED Score (Quantitative Estimate of Drug-likeness)")
        st.caption("QED ranges from 0 (least drug-like) to 1 (most drug-like). Typical approved drugs score 0.5–0.9.")

        fig_gauge, ax_gauge = plt.subplots(figsize=(5, 2.8))
        fig_gauge.patch.set_facecolor("#0e1117")
        ax_gauge.set_facecolor("#0e1117")

        # Background bar
        ax_gauge.barh(0, 1, color="#1e2130", height=0.5, left=0)
        # Color segments
        for start, end, color in [(0, 0.3, "#ef4444"), (0.3, 0.6, "#f59e0b"), (0.6, 1.0, "#10b981")]:
            ax_gauge.barh(0, end - start, color=color, height=0.5, left=start, alpha=0.35)
        # Value bar
        bar_color = "#ef4444" if qed < 0.3 else "#f59e0b" if qed < 0.6 else "#10b981"
        ax_gauge.barh(0, qed, color=bar_color, height=0.5, left=0)
        # Marker line
        ax_gauge.axvline(qed, color="white", linewidth=2.5, ymin=0.15, ymax=0.85)
        # Labels
        ax_gauge.text(qed, 0.38, f"  {qed:.3f}", color="white", fontsize=13, fontweight="bold", va="bottom")
        for x, label in [(0.15, "Low"), (0.45, "Medium"), (0.8, "High")]:
            ax_gauge.text(x, 0, label, color="white", fontsize=8, ha="center", va="center", alpha=0.7)
        ax_gauge.set_xlim(0, 1)
        ax_gauge.set_ylim(-0.4, 0.6)
        ax_gauge.axis("off")
        st.pyplot(fig_gauge, use_container_width=True)
        plt.close(fig_gauge)

        if qed >= 0.6:
            st.success(f"QED = {qed:.3f} — High drug-likeness")
        elif qed >= 0.3:
            st.warning(f"QED = {qed:.3f} — Moderate drug-likeness")
        else:
            st.error(f"QED = {qed:.3f} — Low drug-likeness")

        st.divider()

        # ── Radar chart ─────────────────────────────────────────────────────
        st.subheader("Property Radar vs Drug-likeness Thresholds")
        st.caption("Values normalized against their upper limits. Inside the green zone = passes.")

        labels   = ["MW\n(≤500)", "LogP\n(≤5)", "HBD\n(≤5)", "HBA\n(≤10)", "TPSA\n(≤140)", "RotB\n(≤10)"]
        limits   = [500, 5, 5, 10, 140, 10]
        actuals  = [mw, logp, hbd, hba, tpsa, rb]
        norm     = [min(a / l, 2.0) for a, l in zip(actuals, limits)]

        N = len(labels)
        angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
        norm_plot = norm + norm[:1]
        angles_plot = angles + angles[:1]
        threshold = [1.0] * N + [1.0]

        fig_radar, ax_radar = plt.subplots(figsize=(5, 5), subplot_kw=dict(polar=True))
        fig_radar.patch.set_facecolor("#0e1117")
        ax_radar.set_facecolor("#0e1117")

        # Threshold circle
        ax_radar.plot(angles_plot, threshold, color="#10b981", linewidth=1.5, linestyle="--", alpha=0.8)
        ax_radar.fill(angles_plot, threshold, alpha=0.12, color="#10b981")

        # Molecule polygon
        mol_color = "#10b981" if all(v <= 1 for v in norm) else "#ef4444"
        ax_radar.plot(angles_plot, norm_plot, color=mol_color, linewidth=2)
        ax_radar.fill(angles_plot, norm_plot, alpha=0.3, color=mol_color)

        # Dots at each vertex
        ax_radar.scatter(angles, norm, color=mol_color, s=50, zorder=5)

        ax_radar.set_xticks(angles)
        ax_radar.set_xticklabels(labels, color="white", size=9)
        ax_radar.set_ylim(0, 2.1)
        ax_radar.set_yticks([0.5, 1.0, 1.5, 2.0])
        ax_radar.set_yticklabels(["50%", "limit", "150%", "200%"], color="#6b7280", size=7)
        ax_radar.tick_params(colors="white")
        ax_radar.spines["polar"].set_color("#2e3347")
        ax_radar.grid(color="#2e3347", linewidth=0.8)

        green_patch = mpatches.Patch(color="#10b981", alpha=0.5, label="Threshold zone")
        mol_patch   = mpatches.Patch(color=mol_color, alpha=0.5, label="This molecule")
        ax_radar.legend(handles=[green_patch, mol_patch], loc="upper right",
                        bbox_to_anchor=(1.35, 1.15), labelcolor="white",
                        facecolor="#1e2130", edgecolor="#2e3347", fontsize=8)

        st.pyplot(fig_radar, use_container_width=True)
        plt.close(fig_radar)

        st.divider()

        # ── Scorecard ────────────────────────────────────────────────────────
        st.subheader("Rules Scorecard")
        all_rules = [
            ("Lipinski", "MW ≤ 500 Da",          mw <= 500,   f"{mw:.1f} Da"),
            ("Lipinski", "LogP ≤ 5",              logp <= 5,   f"{logp:.2f}"),
            ("Lipinski", "H-Bond Donors ≤ 5",     hbd <= 5,    str(hbd)),
            ("Lipinski", "H-Bond Acceptors ≤ 10", hba <= 10,   str(hba)),
            ("Veber",    "TPSA ≤ 140 Å²",         tpsa <= 140, f"{tpsa:.1f} Å²"),
            ("Veber",    "Rotatable Bonds ≤ 10",   rb <= 10,    str(rb)),
        ]
        for ruleset, rule, ok, value in all_rules:
            icon  = "✅" if ok else "❌"
            badge = f"`{ruleset}`"
            st.markdown(f"{icon} **{rule}** — {value} &nbsp; {badge}", unsafe_allow_html=True)

        total_violations = sum(1 for _, _, ok, _ in all_rules if not ok)
        st.divider()
        if total_violations == 0:
            st.success("Passes all Lipinski + Veber rules — strong oral drug candidate.")
        elif total_violations <= 1:
            st.warning(f"{total_violations} violation — borderline drug-like.")
        else:
            st.error(f"{total_violations} violations — poor predicted bioavailability.")
