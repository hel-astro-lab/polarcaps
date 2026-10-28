"""
Figure 5-like plot where each curve is replaced by a band covering a range of
collision angles psi (Compton) = theta (two-photon), e.g. 90-180 deg ("head-on regime").

Usage: python plot_lmfp_bands.py msp 90 180
"""
import sys, runpy
import numpy as np

setup = sys.argv[1] if len(sys.argv) > 1 else "msp"
psi_min = float(sys.argv[2]) if len(sys.argv) > 2 else 90.0
psi_max = float(sys.argv[3]) if len(sys.argv) > 3 else 180.0
n_psi = 25

# run the main script to get all the setup quantities
sys.argv = ["estimate_pulsar_quantities.py", setup, "angle", "180"]
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path("estimate_pulsar_quantities.py")

from estimation_functions import comp_eps_ic, comp_lmfp_1phot, comp_lmfp_2phot_angle

tkev, x2, r_pcap = g["tkev"], g["x2"], g["r_pcap"]
gam_smallest, gam_scaling = g["gam_smallest"], g["gam_scaling"]
linestyles, labels = g["linestyles"], g["labels"]
bratio, rad_curv, period = g["bratio"], g["rad_curv"], g["period"]
keV2Kelvin = g["keV2Kelvin"]

psis = np.radians(np.linspace(psi_min, psi_max, n_psi))

res = []
for gs in gam_scaling:
    l1, l2 = [], []
    for psi in psis:
        x1 = comp_eps_ic(gs*gam_smallest, x2, psi)
        l1.append(comp_lmfp_1phot(x1, bratio, rad_curv, period)[0]/r_pcap)
        l2.append(comp_lmfp_2phot_angle(x1, tkev, psi)/r_pcap)
    res.append((np.array(l1), np.array(l2)))

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import NullFormatter, LogLocator

plt.rc("text", usetex=True)
lbfontsz, lwidth = 36, 2.0
lgfontsz = 32   # legend and label-box font size (panels are shown at ~1/4 size in the paper)
plt.rc("xtick", labelsize=lbfontsz); plt.rc("ytick", labelsize=lbfontsz); plt.rc("axes", linewidth=lwidth)
fig, ax = plt.subplots(figsize=(10, 10), dpi=300)

T = tkev*keV2Kelvin
ymax = 1e5 if setup != "wd" else 1e4
for ig in range(len(gam_scaling)):
    for arr, col in ((res[ig][1], "g"), (res[ig][0], "b")):
        a = np.where(np.isfinite(arr), arr, 1e30)
        lo, hi = a.min(axis=0), a.max(axis=0)
        # shade only where every angle in the range is above threshold (finite mfp);
        # near threshold the band is then bounded by the edge curves themselves
        ok = hi < 1e29
        ax.fill_between(T, lo, hi, where=ok, color=col, alpha=0.18, lw=0)
        # edges: psi_max (solid style of the gamma level), psi_min (thin)
        i_hi = np.argmax(psis)
        ax.loglog(T, a[i_hi], color=col, linestyle=linestyles[ig], linewidth=lwidth)
        i_lo = np.argmin(psis)
        ax.loglog(T, a[i_lo], color=col, linestyle=linestyles[ig], linewidth=0.8)

ax.axhline(1.0, color="k", linewidth=2*lwidth, zorder=1)  # l_mfp = H_gap reference
ax.set_xlabel(r"$T_{\mathrm{BB}} [\mathrm{K}]$", fontsize=lbfontsz)
ax.set_ylabel(r"$l_{\mathrm{mfp}}/H_{\mathrm{gap}}$", fontsize=lbfontsz)
if setup == "wd":
    ax.set_xlim(5e4, 1e6); ax.set_ylim(1e-2, 1e4)
else:
    ax.set_xlim(8e4, 1e7); ax.set_ylim(1e-3, 1e5)
secax = ax.secondary_xaxis("top", functions=(lambda K: K/keV2Kelvin, lambda k: k*keV2Kelvin))
secax.set_xlabel(r"$k_\mathrm{B}T_{\mathrm{BB}} [\mathrm{keV}]$", fontsize=lbfontsz, labelpad=12)
if setup == "wd":
    secax.set_xticks([0.005, 0.01, 0.02, 0.05]); secax.set_xticklabels(["$0.005$", "$0.01$", "$0.02$", "$0.05$"])
else:
    secax.set_xticks([0.01, 0.03, 0.1, 0.2, 0.5]); secax.set_xticklabels(["$0.01$", "$0.03$", "$0.1$", "$0.2$", "$0.5$"])
ax.yaxis.set_major_locator(LogLocator(base=10, numticks=20))
ax.xaxis.set_minor_formatter(NullFormatter()); secax.xaxis.set_minor_formatter(NullFormatter())
ax.tick_params(axis="both", which="both", width=lwidth, length=5.0)
secax.tick_params(axis="both", which="both", width=lwidth, length=5.0)

# line-style legend only in the top row (psi up to 180 deg); the styles are the same in all rows
show_legend = (psi_max == 180)
handles = [Line2D([0], [0], color="k", linestyle=linestyles[i], linewidth=lwidth, label=labels[i]) for i in range(len(labels))]
if show_legend:
    if setup == "wd":   # upper right is occupied by the one-photon bands in the WD case
        ax.legend(handles=handles, fontsize=lgfontsz, loc="lower left", handlelength=1.6, labelspacing=0.25, borderpad=0.3)
    else:
        ax.legend(handles=handles, fontsize=lgfontsz, loc="upper right", handlelength=1.6, labelspacing=0.25, borderpad=0.3)
# label box: angle range (thick/thin curves are explained in the caption)
box = ("MSP" if setup != "wd" else "WD") + r", $%g^\circ$--$%g^\circ$" % (psi_min, psi_max)
if setup == "wd" and show_legend:
    bx, by, va = 0.04, 0.96, "top"
else:
    bx, by, va = 0.04, 0.04, "bottom"
ax.text(bx, by, box, transform=ax.transAxes, fontsize=lgfontsz, va=va,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
plt.tight_layout()
fig.savefig("lmfp_band_%s_psi%g-%g.png" % (setup.upper(), psi_min, psi_max), dpi=300)
