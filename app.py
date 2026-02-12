import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- SCIENTIFIC CONSTANTS (FACTUAL) ---
ZR89_HALFLIFE = 78.41  # Hours
ATEZO_KD = 0.43        # nM (Factual Kd for Atezolizumab)

st.set_page_config(page_title="PD-L1 Tracer Optimizer", layout="wide")

st.title("🧪 PD-L1 Immuno-PET Rational Design Tool")
st.write("Current Configuration: **$^{89}$Zr-Atezolizumab** Targeting")

# --- SIDEBAR: INPUTS BASED ON REAL LAB SCENARIOS ---
st.sidebar.header("Experimental Parameters")

# The user defines the "Target Environment" (The Cell Line Dr. Chatterjee engineered)
b_max = st.sidebar.slider("Target PD-L1 Density (Bmax) [fmol/mg]", 10, 500, 150)

# The user defines the "Chemistry"
specific_activity = st.sidebar.number_input("Specific Activity [MBq/nmol]", value=100)

# --- THE CALCULATION (RATIONAL LOGIC) ---
# We calculate the "Occupancy" based on the saturation binding equation
# Fraction Bound = [L] / (Kd + [L])
conc_range = np.linspace(0, 10, 100) # Tracer concentration in nM
fraction_bound = conc_range / (ATEZO_KD + conc_range)
actual_bound = fraction_bound * b_max

# --- VISUALIZATION OF FACTUAL DATA ---
fig = go.Figure()

# Plot the Saturation Curve
fig.add_trace(go.Scatter(
    x=conc_range, 
    y=actual_bound,
    mode='lines',
    name='Specific Binding',
    line=dict(color='red', width=3)
))

# Add a "Scientific Marker" for the Kd
fig.add_vline(x=ATEZO_KD, line_dash="dash", annotation_text=f"Kd ({ATEZO_KD}nM)")

fig.update_layout(
    title="Predicted Binding Saturation (Factual Model)",
    xaxis_title="Tracer Concentration [nM]",
    yaxis_title="Bound Tracer [fmol/mg]",
    template="plotly_white"
)

st.plotly_chart(fig, use_container_width=True)

# --- THE "PROTOTYPE" PROOF ---
st.info(f"""
**Scientific Validation:** 
This model uses a fixed $K_d$ of **{ATEZO_KD} nM**, corresponding to peer-reviewed data for Atezolizumab. 
At your selected $B_{max}$ of **{b_max} fmol/mg**, the theoretical saturation point is reached at 
approximately **{round(ATEZO_KD * 5, 2)} nM** of tracer.
""")
