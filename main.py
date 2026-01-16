import streamlit as st
import os
import sys
import requests
from rdkit import Chem
from rdkit.Chem import Draw

# --- 1. システム・環境設定 ---
# AIエンジンが古いKerasを要求するため設定
os.environ['TF_USE_LEGACY_KERAS'] = '1'

@st.cache_resource
def load_stout_engine():
    try:
        # tensorflow-cpu==2.12.0環境でのロード
        import tensorflow as tf
        from STOUT import translate_forward
        return translate_forward
    except Exception as e:
        return f"AI Load Error: {e}"

# --- 2. 解析・可視化ロジック ---
def show_structure(smiles):
    """SMILESから構造式画像を生成"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        # 旧バージョン互換のため use_column_width を使用
        img = Draw.MolToImage(mol, size=(500, 500))
        st.image(img, use_column_width=True)
        return True
    return False

def fetch_pubchem_all(smiles):
    """PubChemからIUPAC名、慣用名、CIDを取得"""
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
# バージョン固定のため、できるだけシンプルな設定に
st.set_page_config(page_title="StructureEcho", layout="centered")

st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;700&display=swap');
    html, body, [class*="css"] { font-family: 'Inter', sans-serif; }
    .stHeading h1 { font-weight: 800; letter-spacing: -0.05em; margin-bottom: 0px; }
    .result-card { padding: 24px; border-radius: 6px; background-color: #161b22; border: 1px solid #30363d; margin-top: 10px; }
    code { color: #4dabff; background-color: transparent; padding: 0; }
    </style>
    """, unsafe_allow_html=True)

st.title("StructureEcho")
st.caption("MOLECULAR NOMENCLATURE INTELLIGENCE")

stout_model = load_stout_engine()

# クラウド互換性のため radio に変更 (segmented_controlは1.40以降のため)
mode = st.radio(
    "SEARCH ENGINE",
    options=["DATABASE", "AI INFERENCE"],
    horizontal=True
)

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
            # 旧バージョン互換のため use_column_width=True
            if st.button("RUN ANALYSIS", type="primary"):
                with st.spinner("Analyzing..."):
                    if mode == "DATABASE":
                        res = fetch_pubchem_all(smiles_input)
                        if res == "error":
                            st.error("Connection failed.")
                        elif res:
                            st.markdown(f"""
                            <div class="result-card">
                                <p style='color: #8b949e; margin-bottom: 2px; font-size: 0.8rem;'>Common Name</p>
                                <h3 style='margin-top: 0; color: #ffffff;'>{res['common']}</h3>
                                <hr style='border-color: #30363d; margin: 15px 0;'>
                                <p style='color: #8b949e; margin-bottom: 2px; font-size: 0.8rem;'>IUPAC Name</p>
                                <p style='word-break: break-all;'><code>{res['iupac']}</code></p>
                                <div style='margin-top: 15px;'>
                                    <a href="https://pubchem.ncbi.nlm.nih.gov/compound/{res['cid']}" target="_blank" style='font-size: 0.85rem;'>Open Pubchem</a>
                                </div>
                            </div>
                            """, unsafe_allow_html=True)
                        else:
                            st.warning("Not found.")
                    else:
                        if callable(stout_model):
                            try:
                                ai_name = stout_model(smiles_input)
                                st.markdown(f"""
                                <div class="result-card">
                                    <p style='color: #4dabff; margin-bottom: 2px; font-size: 0.8rem;'>AI Predicted Name</p>
                                    <h4 style='margin-top: 0; line-height: 1.4;'>{ai_name}</h4>
                                </div>
                                """, unsafe_allow_html=True)
                            except Exception as e:
                                st.error(f"Inference error: {e}")
                        else:
                            st.error("AI Engine offline.")

with st.sidebar:
    st.markdown("##### System Status")
    if callable(stout_model):
        st.success("STOUT: READY")
    else:
        st.warning("STOUT: LOADING")
    st.divider()
    st.caption("Ver. 3.4.1 (Cloud Optimized)")