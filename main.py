import streamlit as st
import os
import sys
import requests
from rdkit import Chem
from rdkit.Chem import Draw

# --- 1. システム・環境設定 ---
os.environ['TF_USE_LEGACY_KERAS'] = '1'

# 1.11.0環境でのキャッシュ命令
@st.experimental_singleton
def load_stout_engine():
    try:
        import tensorflow as tf
        from STOUT import translate_forward
        return translate_forward
    except Exception as e:
        return f"AI Load Error: {e}"

# --- 2. 解析・可視化ロジック ---
def show_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(500, 500))
        # 1.11.0では use_column_width を使用
        st.image(img, use_column_width=True)
        return True
    return False

def fetch_pubchem_all(smiles):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{requests.utils.quote(smiles)}/property/IUPACName/JSON"
        r = requests.get(url, timeout=5)
        if r.status_code == 200:
            data = r.json()['PropertyTable']['Properties'][0]
            cid = data.get('CID')
            iupac = data.get('IUPACName')
            syn_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
            rs = requests.get(syn_url, timeout=5)
            common = iupac
            if rs.status_code == 200:
                common = rs.json()['InformationList']['Information'][0]['Synonym'][0]
            return {"iupac": iupac, "common": common, "cid": cid}
        return None
    except:
        return "error"

# --- 3. UI 構築 ---
st.set_page_config(page_title="StructureEcho", layout="centered")

# CSSで最低限の整理
st.markdown("""
    <style>
    .result-card { padding: 20px; border-radius: 6px; background-color: #161b22; border: 1px solid #30363d; }
    code { color: #4dabff; }
    </style>
    """, unsafe_allow_html=True)

st.title("StructureEcho")
st.caption("MOLECULAR NOMENCLATURE INTELLIGENCE")

stout_model = load_stout_engine()

# ラジオボタンに戻しました（安定性重視）
mode = st.radio("SEARCH ENGINE", options=["DATABASE", "AI INFERENCE"], horizontal=True)

smiles_input = st.text_input("SMILES String", placeholder="e.g. C1=CC=CC=C1")

if smiles_input:
    col1, col2 = st.columns([1, 1.2])
    with col1:
        st.markdown("##### Structure View")
        valid = show_structure(smiles_input)
    with col2:
        st.markdown("##### Analysis Result")
        if not valid:
            st.error("Invalid SMILES format.")
        else:
            if st.button("RUN ANALYSIS"):
                with st.spinner("Analyzing..."):
                    if mode == "DATABASE":
                        res = fetch_pubchem_all(smiles_input)
                        if res == "error": st.error("Connection failed.")
                        elif res:
                            st.markdown(f"""
                            <div class="result-card">
                                <p style='color: #8b949e; font-size: 0.8rem;'>IUPAC Name</p>
                                <p><code>{res['iupac']}</code></p>
                                <a href="https://pubchem.ncbi.nlm.nih.gov/compound/{res['cid']}" target="_blank">Open Pubchem</a>
                            </div>
                            """, unsafe_allow_html=True)
                        else: st.warning("Not found.")
                    else:
                        if callable(stout_model):
                            try:
                                ai_name = stout_model(smiles_input)
                                st.markdown(f"""
                                <div class="result-card">
                                    <p style='color: #4dabff; font-size: 0.8rem;'>AI Predicted Name</p>
                                    <h4>{ai_name}</h4>
                                </div>
                                """, unsafe_allow_html=True)
                            except Exception as e:
                                st.error(f"Inference error: {e}")
                        else:
                            st.error(f"AI Engine is initializing... Please wait.")

with st.sidebar:
    st.markdown("##### System Status")
    if callable(stout_model):
        st.success("STOUT: READY")
    else:
        st.warning("STOUT: LOADING")
    st.markdown("---")
    st.caption("Ver. 3.4.5 (Stability Focus)")