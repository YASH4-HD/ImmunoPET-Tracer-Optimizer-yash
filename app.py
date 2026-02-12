import streamlit as st
import numpy as np
import plotly.graph_objects as go

# --- SCIENTIFIC CONSTANTS ---
ZR89_HALFLIFE = 78.41  # Hours (Zirconium-89)

TRACER_LIBRARY = {
    "Atezolizumab (Antibody)": {
        "kd_nm": 0.43,
        "molecular_weight_kda": 145,
        "note": "Reference clinical immuno-PET antibody tracer",
    },
    "PD-L1 Nanobody": {
        "kd_nm": 2.10,
        "molecular_weight_kda": 15,
        "note": "Small biologic with fast kinetics and lower MW",
    },
    "PD-L1 Small Molecule": {
        "kd_nm": 12.00,
        "molecular_weight_kda": 0.8,
        "note": "Fast distribution; generally lower affinity than antibodies",
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
**Focus:** PD-L1 Checkpoint Imaging Optimization  
Use mechanistic affinity + occupancy + radioactive decay for protocol tuning.
---
"""
)

# --- SIDEBAR: EXPERIMENTAL DESIGN ---
st.sidebar.header("🧪 Experimental Parameters")

st.sidebar.subheader("Tracer Library")
selected_tracer = st.sidebar.selectbox(
    "Tracer Type",
    options=list(TRACER_LIBRARY.keys()),
    help="Choose tracer class to update affinity (Kd) and molecular weight assumptions.",
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
    help="Target density of PD-L1 receptors in tumor tissue.",
)

uncertainty_pct = st.sidebar.slider(
    "Bmax Uncertainty [%]",
    min_value=0,
    max_value=50,
    value=20,
    help="Creates a confidence band (shaded area) around the binding prediction.",
)

st.sidebar.subheader("Radiochemistry")
spec_act = st.sidebar.number_input(
    "Specific Activity [MBq/nmol]",
    value=100.0,
    min_value=1.0,
    help="Radioactivity per unit of tracer mass.",
)

injected_dose_nm = st.sidebar.slider("Planned Tracer Conc. [nM]", 0.1, 10.0, 2.0)

st.sidebar.subheader("Imaging Logistics")
time_post_injection_h = st.sidebar.slider(
    "Time Post-Injection [hours]",
    min_value=0,
    max_value=120,
    value=48,
    help="Radioactive signal decays over time using Zr-89 half-life.",
)

# --- MODEL (MATH) ---
conc_range = np.linspace(0, 10, 200)

# Central estimate
specific_binding = (b_max * conc_range) / (active_kd + conc_range)
current_binding = (b_max * injected_dose_nm) / (active_kd + injected_dose_nm)

# Uncertainty band from Bmax variance
bmax_low = b_max * (1 - uncertainty_pct / 100)
bmax_high = b_max * (1 + uncertainty_pct / 100)
specific_binding_low = (bmax_low * conc_range) / (active_kd + conc_range)
specific_binding_high = (bmax_high * conc_range) / (active_kd + conc_range)

# Signal conversion and decay
initial_signal_mbq_mg = (current_binding / 1_000_000) * spec_act
decay_constant = np.log(2) / ZR89_HALFLIFE
remaining_fraction = np.exp(-decay_constant * time_post_injection_h)
decayed_signal_mbq_mg = initial_signal_mbq_mg * remaining_fraction

occupancy = (injected_dose_nm / (active_kd + injected_dose_nm)) * 100

# --- VISUALIZATION ---
col1, col2 = st.columns([3, 1])

with col1:
    fig = go.Figure()

    # Lower bound first
    fig.add_trace(
        go.Scatter(
            x=conc_range,
            y=specific_binding_low,
            mode="lines",
            line=dict(width=0),
            showlegend=False,
            hoverinfo="skip",
            name="Lower Uncertainty Bound",
        )
    )

    # Upper bound with fill to previous trace
    fig.add_trace(
        go.Scatter(
            x=conc_range,
            y=specific_binding_high,
            mode="lines",
            line=dict(width=0),
            fill="tonexty",
            fillcolor="rgba(255, 75, 75, 0.20)",
            name=f"Uncertainty Band (±{uncertainty_pct}% Bmax)",
        )
    )

    # Central estimate
    fig.add_trace(
        go.Scatter(
            x=conc_range,
            y=specific_binding,
            mode="lines",
            name="Specific Binding (central estimate)",
            line=dict(color="#FF4B4B", width=4),
        )
    )

    # Operating point
    fig.add_trace(
        go.Scatter(
            x=[injected_dose_nm],
            y=[current_binding],
            mode="markers",
            name="Planned Dose",
            marker=dict(color="black", size=12, symbol="cross"),
        )
    )

    # Kd reference
    fig.add_vline(
        x=active_kd,
        line_dash="dash",
        line_color="green",
        annotation_text=f"Kd ({active_kd} nM)",
    )

    fig.update_layout(
        title="<b>Tracer Binding Saturation Profile with Uncertainty</b>",
        xaxis_title="Tracer Concentration [nM]",
        yaxis_title="Bound Fraction [fmol/mg]",
        hovermode="x unified",
        template="plotly_white",
        height=520,
    )
    st.plotly_chart(fig, width="stretch")

with col2:
    st.metric("Tracer", selected_tracer)
    st.metric("Predicted Binding", f"{current_binding:.2f} fmol/mg")
    st.metric("Signal @ Injection", f"{initial_signal_mbq_mg:.6f} MBq/mg")
    st.metric(
        f"Signal @ {time_post_injection_h}h",
        f"{decayed_signal_mbq_mg:.6f} MBq/mg",
        delta=f"{(remaining_fraction * 100):.1f}% remaining",
    )

    st.progress(occupancy / 100)
    st.write(f"**Receptor Occupancy:** {occupancy:.1f}%")

# --- SCIENTIFIC JUSTIFICATION ---
st.markdown("### 🧠 Scientific Rationale")
st.success(
    f"""
**Validation Logic:**
1. **Affinity Matching:** Active tracer selection updates model affinity using tracer-specific $K_d$ values.
2. **Uncertainty Awareness:** A shaded confidence band (±{uncertainty_pct}% on Bmax) captures biological variability.
3. **Time-Resolved Logistics:** Radioactive decay is explicitly modeled using $N(t)=N_0e^{{-\\lambda t}}$ and Zr-89 half-life ({ZR89_HALFLIFE} h).
4. **Translational Utility:** Dose, occupancy, and practical scan timing can be tuned together to find a realistic imaging sweet spot.
"""
)

st.caption("Developed for: Theranostic Imaging & Immuno-Oncology Research Applications")
