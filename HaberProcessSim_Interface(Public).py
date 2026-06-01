import matplotlib.pyplot as plt
import streamlit as st

def simulate_Haber(temperature_c, pressure_atm, feed_ratio_str):
#-----gas calculation section-----

    T = temperature_c + 273.15
    R_gas = .08206 #(L * atm)/(mol * K)
    R_energy = 8.314 #(J/(mol * K)

    ratio_prts = feed_ratio_str.split(':')
    n2_ratio_pt = float(ratio_prts[0])
    h2_ratio_pt = float(ratio_prts[1])
    ratio_whole = h2_ratio_pt + n2_ratio_pt
    n2_percent = n2_ratio_pt / ratio_whole
    h2_percent = h2_ratio_pt / ratio_whole

    Concentration_total = (pressure_atm)/(R_gas * T)
    n2 = Concentration_total * n2_percent
    h2 = Concentration_total * h2_percent
    nh3 = 0
    print (f'Initial feed concentrations calculated --- [N2] = {n2}M, [H2] = {h2}M')

# -----thermo calculation section-----
    A_for = 1.0e13
    A_rev = 1.0e20
    Ea_for = 138000 #j/mol
    Ea_rev = 230400 #j/mol

    k_for = A_for *(2.718281828 ** (-(Ea_for/(R_energy * T))))
    k_rev = A_rev *(2.718281828 ** (-(Ea_rev/(R_energy * T))))
    if T < 525:
        dt = .1
    elif pressure_atm > 400 or temperature_c > 500:
        dt = .0000001
    elif pressure_atm > 50 or T > 650:
        dt = .000001
    else:
        dt = .0001
    total_time = 300

    pre_nh3 = 0.0
    tol_ = .00000001
    min_steps = 500
    time = 0
    time_history = []
    n2_history = []
    h2_history = []
    nh3_history = []

    while time <= total_time:
        time_history.append(time)
        n2_history.append(n2)
        h2_history.append(h2)
        nh3_history.append(nh3)

        rate_for = k_for * ((max(0.0 , h2)) ** 3) * (max(0.0 , n2))
        rate_rev = k_rev * (max(0.0 , nh3)) ** 2
        rate_net = rate_for - rate_rev

        n2 -= 1 * rate_net * dt
        h2 -= 3 * rate_net * dt
        nh3 += 2 * rate_net * dt

        current_step = len(time_history)
        if current_step > min_steps:
            nh3_change = abs(nh3 - prev_nh3)
            if nh3_change < tol_:
                print(f'Equilibrium reached at {time:.2} seconds.')
                break
        prev_nh3 = nh3
        time += dt

    return time_history, n2_history, h2_history, nh3_history

def plot_Kinetics(time_history, n2_history, h2_history, nh3_history, temperature_c, pressure_atm, feed_ratio_str):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time_history , n2_history , label='Nitrogen ($N_2$)')
    ax.plot(time_history , h2_history , label='Hydrogen ($H_2$)')
    ax.plot(time_history , nh3_history , label='Ammonia ($NH_3$)')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_title(f'Haber Process at {pressure_atm} atm and {temperature_c} C with {feed_ratio_str} feed ratio')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Concentration (M)')
    ax.grid(True)
    st.pyplot(fig)

st.title('Haber-Bosch Process Simulator')
sim_T = st.sidebar.slider('Temperature (°C)', 150, 550, 450)
sim_P = st.sidebar.slider('Pressure (atm)', 100, 500, 200)
st.sidebar.write("### Feed Ratio (N2 : H2)")
col1, col2 = st.sidebar.columns(2)
with col1:
    n2_val = st.selectbox("N2 Part", options=[1, 2, 3, 4], index=0)
with col2:
    h2_val = st.selectbox("H2 Part", options=[1, 2, 3, 4], index=2)

combined_feed_str = f"{n2_val}:{h2_val}"
t_log, n2_log, h2_log, nh3_log = simulate_Haber(sim_T, sim_P, combined_feed_str)
plot_Kinetics(t_log, n2_log, h2_log, nh3_log, sim_T, sim_P, combined_feed_str)