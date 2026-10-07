import itertools
import numpy as np
import matplotlib.pyplot as plt


def plot_all():

    def finish(fig, name):
        fig.tight_layout()
        fig.savefig(name, dpi=200)

    # 7.1.2 Task 2: y = ce^(λx) for c < 0, c = 0, and c > 0
    x = np.linspace(0, 5, 400)
    labc = {"pos": "c > 0", "zero": "c = 0", "neg": "c < 0"}
    labl = {"pos": "λ > 0", "zero": "λ = 0", "neg": "λ < 0"}

    c_full = {"pos": 0.5 * np.arange(1, 7), "zero": [0.0], "neg": -0.5 * np.arange(1, 7)}
    c_short = {"pos": [0.5, 1.5, 2.5], "zero": [0.0], "neg": [-0.5, -1.5, -2.5]}
    l_full = {"pos": 0.1 * np.arange(1, 7), "zero": [0.0], "neg": -0.1 * np.arange(1, 7)}
    l_short = {"pos": [0.2, 0.5], "zero": [0.0], "neg": [-0.2, -0.5]}

    fig, axes = plt.subplots(3, 3, figsize=(15, 11))
    for i, cs in enumerate(["pos", "zero", "neg"]):
        for j, ls in enumerate(["pos", "zero", "neg"]):
            ax = axes[i, j]
            if cs == "zero" and ls == "zero":
                pairs = [(0.0, 0.0)]
            elif cs == "zero":
                pairs = [(0.0, l) for l in l_full[ls]]
            elif ls == "zero":
                pairs = [(c, 0.0) for c in c_full[cs]]
            else:
                pairs = list(itertools.product(c_short[cs], l_short[ls]))
            for c, l in pairs:
                ax.plot(x, c * np.exp(l * x), label=f"c={c:g}, λ={l:g}")
            if cs == "zero":
                ax.set_ylim(-1, 1)
            ax.axhline(0, color="k", lw=0.5)
            ax.set_title(f"{labc[cs]},  {labl[ls]}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=6, ncol=2)
    fig.suptitle("7.1.2: y(x) = c·e^(λx) for all sign combinations of c and λ")
    finish(fig, "fig_7_1_2_c_lambda_grid.png")

    # Eq. (26): 2y' + y = 0, y(x0) = 1  ->  λ = -1/2, c = e^(x0/2)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for x0 in [0, 1, 2]:
        xs = np.linspace(x0, 8, 300)
        ax.plot(xs, np.exp(-(xs - x0) / 2), label=f"x0 = {x0}  (c = e^({x0}/2) = {np.exp(x0/2):.2f})")
        ax.plot(x0, 1, "ko")
    ax.set_title("Eq. (26): 2y' + y = 0, y(x0)=1  →  y = c·e^(−x/2)")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.grid(alpha=0.3); ax.legend()
    finish(fig, "fig_7_1_2_eq26.png")

    # 7.2.2 Task 2c: mass-spring-damper, x0 = 1, xdot0 = 1
    def mass_spring(m, c, k, x0, v0, t):
        disc = c**2 - 4 * m * k
        if disc > 0:
            l1, l2 = np.roots([m, c, k])
            c1 = (v0 - l2 * x0) / (l1 - l2)
            c2 = x0 - c1
            return c1 * np.exp(l1 * t) + c2 * np.exp(l2 * t)
        if disc == 0:
            l = -c / (2 * m)
            c1, c2 = x0, v0 - l * x0
            return (c1 + c2 * t) * np.exp(l * t)
        a = -c / (2 * m)
        b = np.sqrt(4 * m * k - c**2) / (2 * m)
        A, B = x0, (v0 - a * x0) / b
        return np.exp(a * t) * (A * np.cos(b * t) + B * np.sin(b * t))

    t = np.linspace(0, 10, 600)
    cases = [((1, 3, 1), "(1,3,1): distinct real roots"),
             ((1, 2, 2), "(1,2,2): complex roots"),
             ((1, 2, 1), "(1,2,1): repeated root")]
    fig, ax = plt.subplots(figsize=(8, 5))
    for (m, c, k), name in cases:
        ax.plot(t, mass_spring(m, c, k, 1, 1, t), label=f"(m,c,k) = {name}")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_title("7.2.2: mass-spring-damper, x(0)=1, x'(0)=1")
    ax.set_xlabel("t"); ax.set_ylabel("x(t)"); ax.grid(alpha=0.3); ax.legend()
    finish(fig, "fig_7_2_2_mass_spring.png")

    # 7.2.2 Task 3: y''' + 3y'' + 7y' + 5y = 0, y(0)=0, y'(0)=2, y''(0)=3
    # roots: -1, -1 ± 2i
    xx = np.linspace(0, 8, 600)
    y3 = np.exp(-xx) * (7/4 - (7/4) * np.cos(2 * xx) + np.sin(2 * xx))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(xx, y3, label="y(x) = e^(−x)[7/4 − (7/4)cos2x + sin2x]")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_title("7.2.2 Task 3: third-order homogeneous ODE")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.grid(alpha=0.3); ax.legend()
    finish(fig, "fig_7_2_2_third_order.png")

    # 7.3.2 Task 1, Eq. (35): y'' + y = 2e^(-x), y(0)=0, y'(0)=0
    x1 = np.linspace(0, 10, 600)
    yh1 = -np.cos(x1) + np.sin(x1)
    yp1 = np.exp(-x1)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x1, yh1, "--", label="homogeneous  −cos x + sin x")
    ax.plot(x1, yp1, ":", label="particular  e^(−x)")
    ax.plot(x1, yh1 + yp1, "k", lw=2, label="general  −cos x + sin x + e^(−x)")
    ax.set_title("Eq. (35): y'' + y = 2e^(−x)")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.grid(alpha=0.3); ax.legend()
    finish(fig, "fig_7_3_2_eq35.png")

    # 7.3.2 Task 1, Eq. (36): y'' - y = sin x - e^(2x), y(0)=1, y'(0)=-1
    x2 = np.linspace(0, 3, 600)
    yh2 = 0.75 * np.exp(x2) + (7/12) * np.exp(-x2)
    yp2 = -0.5 * np.sin(x2) - (1/3) * np.exp(2 * x2)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x2, yh2, "--", label="homogeneous  ¾eˣ + (7/12)e^(−x)")
    ax.plot(x2, yp2, ":", label="particular  −½sin x − ⅓e^(2x)")
    ax.plot(x2, yh2 + yp2, "k", lw=2, label="general")
    ax.set_title("Eq. (36): y'' − y = sin x − e^(2x)")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.grid(alpha=0.3); ax.legend()
    finish(fig, "fig_7_3_2_eq36.png")

    plt.show()

plot_all()