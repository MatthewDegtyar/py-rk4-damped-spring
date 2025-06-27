import matplotlib.pyplot as plt

# Runge–Kutta method, RK4, for a damped spring

# start with second-order differential equation for damped harmonic oscillator
# m * d²x/dt² + c * dx/dt + k * x = 0

# to apply RK4, convert it to a system of first order ODEs
# Let v = dx/dt (velocity)
# Then:
# dx/dt = v
# dv/dt = - (c/m) * v - (k/m) * x

# express as a vector valued function y(t) = (x(t), v(t))
# and used RK4 to estimate y(t) over time

# TODO: use taylor series to sanity check

def compute_slopes(x_n, v_n, h, c, m, k):
    k_1 = h * v_n
    l_1 = h * ((-c/m)*v_n - (k/m)*x_n)

    x_2_n = x_n + k_1 / 2
    v_2_n = v_n + l_1 / 2
    k_2 = h * v_2_n
    l_2 = h * ((-c/m)*v_2_n - (k/m)*x_2_n)

    x_3_n = x_n + k_2 / 2
    v_3_n = v_n + l_2 / 2
    k_3 = h * v_3_n
    l_3 = h * ((-c/m)*v_3_n - (k/m)*x_3_n)

    x_4_n = x_n + k_3
    v_4_n = v_n + l_3
    k_4 = h * v_4_n
    l_4 = h * ((-c/m)*v_4_n - (k/m)*x_4_n)

    x_next = x_n + (1/6) * (k_1 + 2*k_2 + 2*k_3 + k_4)
    v_next = v_n + (1/6) * (l_1 + 2*l_2 + 2*l_3 + l_4)

    return x_next, v_next


if __name__ == "__main__":
    h = 0.05  # steps
    c = 0.1   # damping coeff
    m = 1.0   # mass
    k = 0.3   # spring constant

    x = 10.0  # init displacement
    v = 0.0  # init velocity
    t = 0.0 # init time

    t_max = 60
    steps = int(t_max / h)

    # storage for plotting
    times = [t]
    positions = [x]
    velocities = [v]

    for _ in range(steps):
        x, v = compute_slopes(x, v, h, c, m, k)
        t += h
        times.append(t)
        positions.append(x)
        velocities.append(v)

    # displacement vs time
    plt.figure(figsize=(10, 5))
    plt.plot(times, positions, label="Displacement x(t)")
    plt.title("Damped Spring Simulation (RK4)")
    plt.xlabel("Time (s)")
    plt.ylabel("Displacement")
    plt.grid(True)
    plt.legend()
    plt.show()
