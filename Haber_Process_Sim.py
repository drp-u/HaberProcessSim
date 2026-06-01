import matplotlib.pyplot as plt

def simulate_Haber():
    print("--Input Reactor Conditions--")
    temperature_c = float(input("Enter the temperature in Celsius: "))
    pressure_atm = float(input("Enter the pressure in atm: "))
    feed_ratio_str = str(input("Enter the feed ratio (N2:H2): "))

#-----gas calculation section-----

    T = temperature_c + 273.15
    R_gas = .08206 #(L * atm)/(mol * K)
    R_energy = 8.1314 #(J/(mol * K)
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

    if pressure_atm > 50 or T > 650:
        dt = .000001
    else:
        dt = .01
    total_time = 1000000

    pre_nh3 = 0.0
    tol_ = .00000001
    min_steps = 200
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

    return time_history, n2_history, h2_history, nh3_history, temperature_c, pressure_atm, feed_ratio_str

def plot_Kinetics(time_history, n2_history, h2_history, nh3_history, temperature_c, pressure_atm, feed_ratio_str):
    plt.plot(time_history , n2_history , label='Nitrogen ($N_2$)')
    plt.plot(time_history , h2_history , label='Hydrogen ($H_2$)')
    plt.plot(time_history , nh3_history , label='Ammonia ($NH_3$)')
    plt.legend(loc='upper right', fontsize=10)
    plt.title(f'Haber Process at {pressure_atm} atm and {temperature_c} C with {feed_ratio_str} feed ratio')
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration (M)')
    plt.grid(True)
    plt.show()

t_log, n2_log, h2_log, nh3_log, sim_T, sim_P, feed_R = simulate_Haber()
plot_Kinetics(t_log, n2_log, h2_log, nh3_log, sim_T, sim_P, feed_R)