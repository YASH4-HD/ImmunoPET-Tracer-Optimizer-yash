import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- SCIENTIFIC CONSTANTS ---
ZR89_HALFLIFE = 78.41  # Hours (Zirconium-89)
ATEZO_KD = 0.43        # nM (Reference Antibody Blocker)
WL12_KD = 25.0         # nM (Chatterjee Peptide Kd)

TRACER_LIBRARY = {
    "WL12 (Chatterjee Peptide)": {
        "kd_nm": 25.0,
        "molecular_weight_kda": 1.5,
        "note": "Optimized peptide for PD-L1; fast clearance.",
    },
    "Atezolizumab (Antibody)": {
        "kd_nm": 0.43,
        "molecular_weight_kda": 145,
        "note": "Reference clinical immuno-PET antibody tracer.",
    },
    "PD-L1 Nanobody": {
        "kd_nm": 2.10,
        "molecular_weight_kda": 15,
        "note": "Small biologic with fast kinetics.",
    },
}

st.set_page_config(
    page_title="ImmunoPET-Tracer-Optimizer",
    page_icon="🧬",
    layout="wide",
)

# --- HEADER ---
st.title("🎯 Immuno-PET Tracer Rational Design Tool")
st.markdown(
    """
**Advanced Predictive Modeling:** PD-L1 Checkpoint Imaging & Target Engagement  
*Integrating Mechanistic Binding, Competitive Inhibition, and Pharmacokinetic Sink Effects.*
---
"""
)

# --- SIDEBAR: EXPERIMENTAL DESIGN ---
st.sidebar.header("🧪 Experimental Parameters")

st.sidebar.subheader("Tracer Library")
selected_tracer = st.sidebar.selectbox(
    "Select Tracer for Simulation",
    options=list(TRACER_LIBRARY.keys()),
    help="Updates affinity (Kd) and MW assumptions based on tracer class.",
)

active_tracer = TRACER_LIBRARY[selected_tracer]
active_kd = active_tracer["kd_nm"]
active_mw = active_tracer["molecular_weight_kda"]

st.sidebar.info(
    f"**Selected Kd:** {active_kd} nM\n\n"
    f"**Molecular Weight:** {active_mw} kDa\n\n"
    f"{active_tracer['note']}"
)

st.sidebar.subheader("Biological Target")
b_max = st.sidebar.slider(
    "Receptor Density (Bmax) [fmol/mg]",
    min_value=10,
    max_value=500,
    value=150,
)

uncertainty_pct = st.sidebar.slider(
    "Bmax Uncertainty [%]",
    min_value=0,
    max_value=50,
    value=20,
    help="Confidence band for biological variability.",
)

# --- CHATTERJEE SPECIFIC INPUTS ---
st.sidebar.subheader("Competitive Environment")
ab_blocker_conc = st.sidebar.slider(
    "Therapeutic Antibody Blocker [nM]",
    0.0, 20.0, 0.0,
    help="Simulates a patient already on antibody therapy (e.g., Atezolizumab)."
)

st.sidebar.subheader("Pharmacokinetics (The Sink)")
liver_sink = st.sidebar.slider(
    "Liver/Kidney Sequestration [%]",
    0, 90, 20,
    help="Percentage of dose lost to off-target organs (The 'Sink Effect')."
)

st.sidebar.subheader("Radiochemistry & Timing")
spec_act = st.sidebar.number_input("Specific Activity [MBq/nmol]", value=100.0)
injected_dose_nm = st.sidebar.slider("Planned Tracer Injection [nM]", 0.1, 50.0, 5.0)
time_h = st.sidebar.slider("Time Post-Injection [hours]", 0, 120, 24)

# --- MODEL MATH (INTEGRATED LOGIC) ---

# 1. Apply Sink Effect to injected dose
effective_conc = injected_dose_nm * (1 - (liver_sink / 100))
conc_range = np.linspace(0.1, 50, 200)
effective_range = conc_range * (1 - (liver_sink / 100))

# 2. Competitive Binding Equation (The core scientific novelty)
# Kd_app = Kd * (1 + [Blocker]/Kd_blocker)
kd_apparent = active_kd * (1 + (ab_blocker_conc / ATEZO_KD))

def calc_binding(c, bmax_val, kd_app):
    return (bmax_val * c) / (kd_app + c)

# Central, High, and Low estimates
specific_binding = calc_binding(effective_range, b_max, kd_apparent)
b_low, b_high = b_max * (1 - uncertainty_pct/100), b_max * (1 + uncertainty_pct/100)
binding_low = calc_binding(effective_range, b_low, kd_apparent)
binding_high = calc_binding(effective_range, b_high, kd_apparent)

# Current point calculation
current_binding = calc_binding(effective_conc, b_max, kd_apparent)
occupancy = (effective_conc / (kd_apparent + effective_conc)) * 100

# 3. Radioactive Decay
decay_constant = np.log(2) / ZR89_HALFLIFE
remaining_frac = np.exp(-decay_constant * time_h)
initial_signal = (current_binding / 1_000_000) * spec_act
decayed_signal = initial_signal * remaining_frac

# --- VISUALIZATION ---
col1, col2 = st.columns([3, 1])

with col1:
    fig = go.Figure()

    # Uncertainty Band
    fig.add_trace(go.Scatter(x=conc_range, y=binding_low, mode="lines", line=dict(width=0), showlegend=False))
    fig.add_trace(go.Scatter(
        x=conc_range, y=binding_high, mode="lines", line=dict(width=0),
        fill="tonexty", fillcolor="rgba(255, 75, 75, 0.2)",
        name=f"Uncertainty (±{uncertainty_pct}% Bmax)"
    ))

    # Binding Curve
    fig.add_trace(go.Scatter(
        x=conc_range, y=specific_binding,
        mode="lines", name=f"{selected_tracer} Binding",
        line=dict(color="#FF4B4B", width=4)
    ))

    # Operating Point
    fig.add_trace(go.Scatter(
        x=[injected_dose_nm], y=[current_binding],
        mode="markers", name="Operating Point",
        marker=dict(color="black", size=12, symbol="cross")
    ))

    fig.update_layout(
        title=f"<b>Binding Profile: {selected_tracer}</b><br><sup>Effect of {liver_sink}% Sink & {ab_blocker_conc}nM Blocker</sup>",
        xaxis_title="Injected Concentration [nM]",
        yaxis_title="Bound Fraction [fmol/mg]",
        template="plotly_white", height=550
    )
    st.plotly_chart(fig, use_container_width=True)

with col2:
    st.metric("Effective Delivery", f"{effective_conc:.2f} nM", delta=f"-{liver_sink}% Sink")
    st.metric("Predicted Binding", f"{current_binding:.2f} fmol/mg")
    st.metric(f"Signal @ {time_h}h", f"{decayed_signal:.6f} MBq/mg")
    
    st.write(f"**Receptor Occupancy:** {occupancy:.1f}%")
    st.progress(min(occupancy/100, 1.0))
    
    if ab_blocker_conc > 0:
        st.warning(f"Kd shifted to {kd_apparent:.2f} nM due to competition.")

# --- SCIENTIFIC JUSTIFICATION ---
st.markdown("### 🧠 Multi-Factor Validation Logic")
st.success(
    f"""
1. **Competitive Engagement:** Uses the Cheng-Prusoff derived model to calculate $K_{{d,app}}$ when therapeutic antibodies are present.
2. **Pharmacokinetic Sink:** Accounts for the liver/kidney sequestration identified in WL12 biodistribution studies.
3. **Uncertainty Propagation:** Shaded regions represent the impact of heterogeneous PD-L1 expression ($B_{{max}}$) on signal reliability.
4. **Isotope Decay:** Real-time signal attenuation based on $^{{89}}$Zr kinetics ($T_{{1/2}} = 78.41$h).
"""
)

st.caption("Developed for: Chatterjee Lab & Translational Immuno-Oncology Research")
