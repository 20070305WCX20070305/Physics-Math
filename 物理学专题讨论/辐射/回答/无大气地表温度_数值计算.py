# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

S0 = 1361.0
SIGMA = 5.670374419e-8
P = 86400.0
OMEGA = 2.0 * np.pi / P
EMISS = 1.0
DZ = 0.004
DEPTH = 1.2
DT = 20.0
DAYS = 25

MATERIALS = {
    "花岗岩": (2.5, 2.3e6),
    "湿土": (1.5, 2.5e6),
    "干沙": (0.3, 1.2e6),
}


def simulate(lat_deg, k, rhoc, albedo=0.30, days=DAYS, dt=DT, dz=DZ, depth=DEPTH):
    n = int(round(depth / dz)) + 1
    z = np.linspace(0.0, depth, n)
    kap = k / rhoc
    phi = np.radians(lat_deg)
    steps_per_day = int(round(P / dt))
    nsteps = days * steps_per_day
    fbar = (1.0 - albedo) * S0 * np.cos(phi) / np.pi
    t_iso = (fbar / (EMISS * SIGMA)) ** 0.25
    temp = np.full(n, t_iso)
    r = kap * dt / dz ** 2
    cap = rhoc * dz / (2.0 * dt)
    rec_t = np.empty(steps_per_day)
    rec_ts = np.empty(steps_per_day)
    prof = {}
    day_balance = []
    e_abs = 0.0
    e_emit = 0.0
    for step in range(nsteps):
        t_now = (step + 1) * dt
        hour = (t_now % P) / P * 2.0 * np.pi
        mu = max(np.cos(phi) * np.cos(hour), 0.0)
        f_solar = (1.0 - albedo) * S0 * mu
        temp_old = temp.copy()
        old0 = temp_old[0]
        star = old0
        for _ in range(3):
            a0 = cap + k / dz + 4.0 * EMISS * SIGMA * star ** 3
            b0 = -k / dz
            rhs0 = cap * old0 + f_solar + 3.0 * EMISS * SIGMA * star ** 4
            ab = np.empty((3, n))
            ab[0, :] = -r
            ab[1, :] = 1.0 + 2.0 * r
            ab[2, :] = -r
            ab[0, 1] = b0
            ab[1, 0] = a0
            ab[1, -1] = 1.0
            ab[2, -2] = -1.0
            rhs = np.empty(n)
            rhs[0] = rhs0
            rhs[1:-1] = temp_old[1:-1]
            rhs[-1] = 0.0
            sol = solve_banded((1, 1), ab, rhs)
            star = sol[0]
        temp = sol
        e_abs += f_solar * dt
        e_emit += EMISS * SIGMA * temp[0] ** 4 * dt
        if step >= nsteps - steps_per_day:
            j = step - (nsteps - steps_per_day)
            rec_t[j] = t_now % P
            rec_ts[j] = temp[0]
            if j == 0:
                prof["正午"] = temp.copy()
            if j == steps_per_day // 2:
                prof["午夜"] = temp.copy()
        if (step + 1) % steps_per_day == 0:
            day = (step + 1) // steps_per_day
            imbalance = (e_abs - e_emit) / P
            day_balance.append(imbalance)
            if day < days - 1:
                temp += imbalance * P / (rhoc * depth)
            e_abs = 0.0
            e_emit = 0.0
    spec = np.fft.rfft(rec_ts)
    amps = 2.0 * np.abs(spec) / steps_per_day
    t_eff = (np.mean(rec_ts ** 4)) ** 0.25
    beta = np.sqrt(k * rhoc * OMEGA)
    dskin = np.sqrt(2.0 * kap / OMEGA)
    yrad = 4.0 * EMISS * SIGMA * t_eff ** 3
    ymag = np.hypot(yrad + beta / np.sqrt(2.0), beta / np.sqrt(2.0))
    f1 = (1.0 - albedo) * S0 * np.cos(phi) / 2.0
    res = dict(
        lat=lat_deg, albedo=albedo, k=k, rhoc=rhoc,
        kappa=kap, skin=dskin, beta=beta, yrad=yrad, ymag=ymag,
        t_iso=t_iso, t_eff=t_eff,
        t_mean=float(rec_ts.mean()), t_min=float(rec_ts.min()),
        t_max=float(rec_ts.max()), span=float(rec_ts.max() - rec_ts.min()),
        amp1=float(amps[1]), amp2=float(amps[2]), amp3=float(amps[3]),
        dt1_ana=float(f1 / ymag),
        balances=day_balance,
        z=z, rec_t=rec_t, rec_ts=rec_ts, prof=prof,
    )
    return res


def day_curve(res):
    n = len(res["rec_ts"])
    half = n // 2
    y = np.roll(res["rec_ts"], -half)
    x = np.arange(n) * (P / n) / 3600.0
    return x, y


CASES = [
    ("EQ granite", 0.0, "花岗岩", 0.30),
    ("EQ wet soil", 0.0, "湿土", 0.30),
    ("EQ dry sand", 0.0, "干沙", 0.30),
    ("30N dry sand", 30.0, "干沙", 0.30),
    ("45N dry sand", 45.0, "干沙", 0.30),
    ("60N dry sand", 60.0, "干沙", 0.30),
    ("45N granite", 45.0, "花岗岩", 0.30),
    ("EQ dry sand a=.12", 0.0, "干沙", 0.12),
]

store = {}
print(f"{'case':>18} {'T_eff':>7} {'T_mean':>7} {'T_min':>7} {'T_max':>7} {'span':>7} {'amp1':>6} {'amp2':>6} {'T1_ana':>7} {'bal':>8}")
for name, lat, mat, alb in CASES:
    k, rhoc = MATERIALS[mat]
    res = simulate(lat, k, rhoc, albedo=alb)
    store[(lat, mat, alb)] = res
    store[name] = res
    print(f"{name:>18} {res['t_eff']:7.1f} {res['t_mean']:7.1f} {res['t_min']:7.1f} {res['t_max']:7.1f} {res['span']:7.1f} "
          f"{res['amp1']:6.1f} {res['amp2']:6.1f} {res['dt1_ana']:7.1f} {res['balances'][-1]:8.2f}")
    print(f"    skin={res['skin']*100:.1f} cm  beta={res['beta']:.1f} W/m2/K  yrad={res['yrad']:.2f} W/m2/K  ytot={res['ymag']:.1f} W/m2/K")

fig, axes = plt.subplots(1, 3, figsize=(17.0, 4.8), constrained_layout=True)

ax = axes[0]
colors = {"干沙": "tab:red", "湿土": "tab:green", "花岗岩": "tab:blue"}
for mat in ["干沙", "湿土", "花岗岩"]:
    res = store[(0.0, mat, 0.30)]
    x, y = day_curve(res)
    ax.plot(x, y, color=colors[mat], lw=1.8, label=mat)
ax.axhline(store[(0.0, "干沙", 0.30)]["t_iso"], color="k", ls="--", lw=1.0)
ax.text(0.3, store[(0.0, "干沙", 0.30)]["t_iso"] + 4, r"$T_{\rm eff}=270\,$K（日平均辐射平衡）", fontsize=9)
ax.axvspan(6, 18, color="gold", alpha=0.12)
ax.set_xlim(0, 24)
ax.set_xticks([0, 6, 12, 18, 24])
ax.set_xlabel("地方时 / h")
ax.set_ylabel("地表温度 / K")
ax.set_title("(a) 赤道不同地表材料")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

ax = axes[1]
for lat, style in [(0.0, "-"), (30.0, "--"), (45.0, "-."), (60.0, ":")]:
    res = store[(lat, "干沙", 0.30)]
    x, y = day_curve(res)
    ax.plot(x, y, color="tab:red", ls=style, lw=1.8, label=f"{int(lat)}°")
ax.axvspan(6, 18, color="gold", alpha=0.12)
ax.set_xlim(0, 24)
ax.set_xticks([0, 6, 12, 18, 24])
ax.set_xlabel("地方时 / h")
ax.set_ylabel("地表温度 / K")
ax.set_title("(b) 干沙地表随纬度的变化")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

ax = axes[2]
res = store[(0.0, "干沙", 0.30)]
zmask = res["z"] <= 0.8
ax.plot(res["prof"]["正午"][zmask], res["z"][zmask], color="tab:orange", lw=1.8, label="正午")
ax.plot(res["prof"]["午夜"][zmask], res["z"][zmask], color="tab:purple", lw=1.8, label="午夜")
d = res["skin"]
ax.axhline(d, color="gray", ls="--", lw=1.0)
ax.text(d + 0.01, 0.5, f"日周期穿透深度 $d={d*100:.0f}$ cm", fontsize=9, rotation=90, va="center")
ax.invert_yaxis()
ax.set_xlabel("温度 / K")
ax.set_ylabel("深度 z / m")
ax.set_title("(c) 地下温度剖面（赤道干沙）")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

fig.savefig("fig_无大气地表温度.png", dpi=150)
print("figure saved")
