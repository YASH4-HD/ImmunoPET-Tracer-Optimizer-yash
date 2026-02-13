import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- SCIENTIFIC CONSTANTS (UPDATED PER LESNIAK ET AL. 2019) ---
ZR89_HALFLIFE = 78.41  # Hours (Zirconium-89)
ATEZO_KD = 0.43        # nM (Reference Antibody Blocker from 2016 paper)
WL12_KI = 12.3         # nM (Ki from 2019 Paper - more accurate than IC50)

TRACER_LIBRARY = {
    "WL12 (Chatterjee Peptide)": {
        "kd_nm": 12.3,
        "molecular_weight_kda": 1.5,
        "note": "Affinity based on Ki (12.3nM) from 2019 FRET assay. Accounts for FPy-WL12 analog.",
    },
    "Atezolizumab (Antibody)": {
        "kd_nm": 0.43,
        "molecular_weight_kda": 145,
        "note": "Reference clinical immuno-PET antibody tracer (Sub-nanomolar affinity).",
    },
    "PD-L1 Nanobody": {
        "kd_nm": 2.10,
        "molecular_weight_kda": 15,
        "note": "Small biologic with fast kinetics and mid-range affinity.",
    },
    "PD-L1 Small Molecule": {
        "kd_nm": 12.00,
        "molecular_weight_kda": 0.8,
        "note": "Fast distribution; generally lower affinity than antibodies.",
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
*Integrating Mechanistic Binding (Cheng-Prusoff), Pharmacokinetic Sink Effects, and Isotope Decay.*
---
"""
)

# --- SIDEBAR: EXPERIMENTAL DESIGN ---
st.sidebar.header("🧪 Experimental Parameters")

st.sidebar.subheader("Tracer Library")
selected_tracer = st.sidebar.selectbox(
    "Select Tracer for Simulation",
    options=list(TRACER_LIBRARY.keys()),
    help="Updates affinity (Kd/Ki) and MW assumptions based on tracer class.",
)

active_tracer = TRACER_LIBRARY[selected_tracer]
active_kd = active_tracer["kd_nm"]
active_mw = active_tracer["molecular_weight_kda"]

st.sidebar.info(
    f"**Selected Kd/Ki:** {active_kd} nM\n\n"
    f"**Molecular Weight:** {active_mw} kDa\n\n"
    f"{active_tracer['note']}"
)

st.sidebar.subheader("Biological Target")
b_max = st.sidebar.slider(
    "Receptor Density (Bmax) [fmol/mg]",
    min_value=10,
    max_value=500,
    value=150,
    help="Target density of PD-L1 receptors in tumor tissue."
)

uncertainty_pct = st.sidebar.slider(
    "Bmax Uncertainty [%]",
    min_value=0,
    max_value=50,
    value=20,
    help="Creates a confidence band (shaded area) around the binding prediction.",
)

# --- COMPETITION & SINK (THE CHATTERJEE VARIABLES) ---
st.sidebar.subheader("Competitive Environment")
ab_blocker_conc = st.sidebar.slider(
    "Therapeutic Antibody Blocker [nM]",
    0.0, 20.0, 0.0,
    help="Simulates a patient already on antibody therapy (e.g., Atezolizumab)."
)

st.sidebar.subheader("Pharmacokinetics (The Sink)")
liver_sink = st.sidebar.slider(
    "Liver/Kidney Sequestration [%]",
    0, 90, 30,
    help="Total % of injected dose lost to off-target organs. (Based on ~20% ID/g in 2019 Paper)."
)
st.sidebar.caption("Note: Sink % represents total dose unavailable for tumor binding.")

st.sidebar.subheader("Radiochemistry & Timing")
spec_act = st.sidebar.number_input(
    "Specific Activity [MBq/nmol]", 
    value=100.0,
    help="Radioactivity per unit of tracer mass."
)
injected_dose_nm = st.sidebar.slider("Planned Tracer Injection [nM]", 0.1, 50.0, 5.0)
time_h = st.sidebar.slider(
    "Time Post-Injection [hours]", 
    0, 120, 24,
    help="Radioactive signal decays over time using Zr-89 half-life."
)

# --- MODEL MATH (INTEGRATED LOGIC) ---

# 1. Apply Sink Effect to injected dose
effective_conc = injected_dose_nm * (1 - (liver_sink / 100))
conc_range = np.linspace(0.1, 50, 200)
effective_range = conc_range * (1 - (liver_sink / 100))

# 2. Competitive Binding Equation (Cheng-Prusoff)
# Kd_apparent = Kd * (1 + [Blocker]/Kd_blocker)
kd_apparent = active_kd * (1 + (ab_blocker_conc / ATEZO_KD))

def calc_binding(c, bmax_val, kd_app):
    return (bmax_val * c) / (kd_app + c)

# Central, High, and Low estimates for Uncertainty Band
specific_binding = calc_binding(effective_range, b_max, kd_apparent)
b_low = b_max * (1 - uncertainty_pct / 100)
b_high = b_max * (1 + uncertainty_pct / 100)
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

    # Uncertainty Band (Fill)
    fig.add_trace(go.Scatter(
        x=conc_range, y=binding_low, 
        mode="lines", line=dict(width=0), 
        showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=conc_range, y=binding_high, 
        mode="lines", line=dict(width=0),
        fill="tonexty", fillcolor="rgba(255, 75, 75, 0.2)",
        name=f"Uncertainty (±{uncertainty_pct}% Bmax)"
    ))

    # Central Binding Curve
    fig.add_trace(go.Scatter(
        x=conc_range, y=specific_binding,
        mode="lines", name=f"{selected_tracer} Binding",
        line=dict(color="#FF4B4B", width=4)
    ))

    # Operating Point
    fig.add_trace(go.Scatter(
        x=[injected_dose_nm], y=[current_binding],
        mode="markers", name="Planned Dose",
        marker=dict(color="black", size=12, symbol="cross")
    ))

    # Kd Apparent Marker
    fig.add_vline(
        x=kd_apparent, 
        line_dash="dash", 
        line_color="green",
        annotation_text=f"Kd_app ({kd_apparent:.1f} nM)"
    )

    fig.update_layout(
        title=f"<b>Tracer Binding Profile: {selected_tracer}</b><br><sup>Modeling {liver_sink}% Sink Effect & {ab_blocker_conc}nM Antibody Competition</sup>",
        xaxis_title="Injected Concentration [nM]",
        yaxis_title="Bound Fraction [fmol/mg]",
        hovermode="x unified",
        template="plotly_white", 
        height=550
    )
    st.plotly_chart(fig, use_container_width=True)

with col2:
    st.metric("Effective Delivery", f"{effective_conc:.2f} nM", delta=f"-{liver_sink}% Sink")
    st.metric("Predicted Binding", f"{current_binding:.2f} fmol/mg")
    st.metric("Signal @ Injection", f"{initial_signal:.6f} MBq/mg")
    st.metric(
        f"Signal @ {time_h}h", 
        f"{decayed_signal:.6f} MBq/mg",
        delta=f"{(remaining_frac * 100):.1f}% remaining"
    )
    
    st.write(f"**Receptor Occupancy:** {occupancy:.1f}%")
    st.progress(min(occupancy/100, 1.0))
    
    if ab_blocker_conc > 0:
        st.warning(f"Affinity shifted from {active_kd} to {kd_apparent:.2f} nM due to competition.")

# --- SCIENTIFIC JUSTIFICATION ---
st.markdown("### 🧠 Multi-Factor Validation Logic")
st.success(
    f"""
**Validation Framework:**
1. **Competitive Engagement:** Uses the Cheng-Prusoff derived model to calculate $K_{{d,app}}$ ($K_{{d,app}} = K_i \cdot (1 + [Drug]/K_{{d,drug}})$) when therapeutic antibodies are present.
2. **Pharmacokinetic Sink:** Accounts for the significant hepatic/renal sequestration (~20% ID/g) identified in the 2019 WL12 biodistribution studies.
3. **Uncertainty Propagation:** Shaded regions represent the impact of heterogeneous PD-L1 expression ($B_{{max}}$) on signal reliability.
4. **Translational Logistics:** Real-time signal attenuation based on $^{{89}}$Zr kinetics ($T_{{1/2}} = 78.41$h) allows for optimization of delayed imaging windows.
"""
)

st.caption("Developed for: Chatterjee Lab | Translational Immuno-Oncology Research Applications")
