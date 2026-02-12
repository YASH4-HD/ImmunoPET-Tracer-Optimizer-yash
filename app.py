import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- SCIENTIFIC CONSTANTS (FACTUAL & PEER-REVIEWED) ---
# Based on Dr. Chatterjee's field of Immuno-PET
ZR89_HALFLIFE = 78.41  # Hours (Zirconium-89)
ATEZO_KD = 0.43        # nM (Dissociation constant for Anti-PD-L1 Atezolizumab)

st.set_page_config(
    page_title="ImmunoPET-Tracer-Optimizer",
    page_icon="🧬",
    layout="wide"
)

# --- HEADER ---
st.title("🎯 Immuno-PET Tracer Rational Design Tool")
st.markdown(f"""
**Focus:** PD-L1 Checkpoint Imaging Optimization  
**Reference Tracer:** $^{{89}}$Zr-DFO-Atezolizumab ($K_d \approx {ATEZO_KD}$ nM)
---
""")

# --- SIDEBAR: EXPERIMENTAL DESIGN ---
st.sidebar.header("🧪 Experimental Parameters")

st.sidebar.subheader("Biological Target")
# Bmax range based on recombinant PD-L1 overexpressing cell lines
b_max = st.sidebar.slider(
    "Receptor Density (Bmax) [fmol/mg]", 
    min_value=10, 
    max_value=500, 
    value=150,
    help="Target density of PD-L1 receptors in the tumor tissue."
)

st.sidebar.subheader("Radiochemistry")
spec_act = st.sidebar.number_input(
    "Specific Activity [MBq/nmol]", 
    value=100,
    help="Radioactivity per unit of antibody mass."
)

injected_dose_nm = st.sidebar.slider(
    "Planned Tracer Conc. [nM]", 
    0.1, 10.0, 2.0
)

# --- THE RATIONAL MODEL (MATH) ---
# 1. Saturation Binding Equation (Langmuir Isotherm)
conc_range = np.linspace(0, 10, 200) 
specific_binding = (b_max * conc_range) / (ATEZO_KD + conc_range)

# 2. Predicted Signal at specific injected dose
current_binding = (b_max * injected_dose_nm) / (ATEZO_KD + injected_dose_nm)
# Convert fmol/mg to MBq/mg: (fmol -> nmol = /1e6) * MBq/nmol
predicted_signal = (current_binding / 1000000) * spec_act

# --- VISUALIZATION ---
col1, col2 = st.columns([3, 1])

with col1:
    fig = go.Figure()

    # The Saturation Curve
    fig.add_trace(go.Scatter(
        x=conc_range, 
        y=specific_binding,
        mode='lines',
        name='Specific Binding (fmol/mg)',
        line=dict(color='#FF4B4B', width=4)
    ))

    # Current Operating Point
    fig.add_trace(go.Scatter(
        x=[injected_dose_nm], 
        y=[current_binding],
        mode='markers',
        name='Planned Dose',
        marker=dict(color='black', size=12, symbol='cross')
    ))

    # Kd Reference Line
    fig.add_vline(x=ATEZO_KD, line_dash="dash", line_color="green", 
                  annotation_text=f"Kd ({ATEZO_KD}nM)")

    fig.update_layout(
        title="<b>Tracer Binding Saturation Profile</b>",
        xaxis_title="Tracer Concentration [nM]",
        yaxis_title="Bound Fraction [fmol/mg]",
        hovermode="x unified",
        template="plotly_white",
        height=500
    )
    st.plotly_chart(fig, use_container_width=True)

with col2:
    st.metric("Predicted Binding", f"{current_binding:.2f} fmol/mg")
    st.metric("Imaging Signal", f"{predicted_signal:.6f} MBq/mg")
    
    occupancy = (injected_dose_nm / (ATEZO_KD + injected_dose_nm)) * 100
    st.progress(occupancy/100)
    st.write(f"**Receptor Occupancy:** {occupancy:.1f}%")

# --- SCIENTIFIC JUSTIFICATION ---
st.markdown("### 🧠 Scientific Rationale")
st.success(f"""
**Validation Logic:**
1. **Affinity Matching:** The model uses a fixed $K_d$ of **{ATEZO_KD} nM**, ensuring calculations remain consistent with Atezolizumab's known binding kinetics.
2. **Signal Prediction:** By integrating Specific Activity (**{spec_act} MBq/nmol**), the tool predicts the absolute radioactive signal. This allows for 'Rational Scaling'—ensuring that the injected dose is sufficient for PET detection while avoiding receptor saturation.
3. **Translational Utility:** This helps in selecting the 'Sweet Spot' for imaging—usually between 50% and 80% occupancy—to maximize signal without masking the biological target.
""")

st.caption("Developed for: Theranostic Imaging & Immuno-Oncology Research Applications")
