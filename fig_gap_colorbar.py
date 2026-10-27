"""
Color bars for the phase-space panels of Figure 2 (made by fig_gap_MSP_WD.py).

The panels show log10 of the 2D histograms h2_ene_e-, h2_ene_e+ and h2_ene_ph written by
QEDToolset (qed_toolset.py), with fixed color limits toolset.pylims (pairs) and
toolset.xylims (photons). Here the histogram values are converted to physical units:

  pairs:   h2 = sum(w) / (dlog10 p * N_box)                over one height bin
  photons: h2 = sum(w eps) / (dlog10 eps * N_time * N_box) over one height bin

where a height bin contains (Nx*NxMesh)/Nhist cells and the co-rotation density n_co corresponds
to conf.ppc (weight-1) particles per cell. Hence

  h2_pairs   = f * (dn/dlog10 p) / n_co,
  h2_photons = f * (eps dn_eps/dlog10 eps) / n_co / N_time,
  f = ppc * (cells per height bin) / N_box          (= 16*80/1024 = 1.25 for msp_2ph.ini)

and the color bars below are labeled in these normalized units. For pairs this is the local
multiplicity per decade, dm/dlog10 p (m = n/n_co); e- + e+ summed matches the Fig. 4 spectrum.

Usage:  python3 fig_gap_colorbar.py --conf msp_2ph.ini
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import pytools
from init_problem import Configuration_Pulsar as Configuration
from qed_toolset import QEDToolset

if __name__ == "__main__":
    args = pytools.parse_args()
    conf = Configuration(args.conf_filename, do_print=False)
    fdir = conf.outdir + "/"
    toolset = QEDToolset(conf)

    # conversion factor from histogram values to (dn/dlog10 p)/n_co
    nx = conf.Nx*conf.NxMesh
    cells_per_bin = nx/toolset.Nhist
    f = conf.ppc*cells_per_bin/toolset.N_box
    f_ph = f/toolset.N_time
    print("cells per height bin:", cells_per_bin, " f (pairs):", f, " f (photons):", f_ph)

    plt.rc('font', family='serif')
    plt.rc('text', usetex=True)
    plt.rc('xtick', labelsize=7)
    plt.rc('axes', labelsize=8)

    fig = plt.figure(figsize=(7.0, 0.75))
    ax1 = fig.add_axes([0.08, 0.55, 0.38, 0.18])   # pairs
    ax2 = fig.add_axes([0.56, 0.55, 0.38, 0.18])   # photons

    cmap = matplotlib.colormaps['turbo']

    # same color limits as in fig_gap_MSP_WD.py, shifted to normalized units
    n1 = Normalize(vmin=np.log10(toolset.pylims[0]/f),    vmax=np.log10(toolset.pylims[1]/f))
    n2 = Normalize(vmin=np.log10(toolset.xylims[0]/f_ph), vmax=np.log10(toolset.xylims[1]/f_ph))

    cb1 = fig.colorbar(matplotlib.cm.ScalarMappable(norm=n1, cmap=cmap), cax=ax1, orientation='horizontal')
    cb2 = fig.colorbar(matplotlib.cm.ScalarMappable(norm=n2, cmap=cmap), cax=ax2, orientation='horizontal')

    cb1.set_ticks(np.arange(np.ceil(n1.vmin), n1.vmax, 1.0))
    cb2.set_ticks(np.arange(np.ceil(n2.vmin), n2.vmax, 1.0))
    cb1.minorticks_on(); cb2.minorticks_on()

    cb1.set_label(r"$\log_{10}\left(\mathrm{d}m/\mathrm{d}\log_{10} p\right)$ \ (pairs)")
    cb2.set_label(r"$\log_{10}\left[(\varepsilon\,\mathrm{d}n_\varepsilon/\mathrm{d}\log_{10} \varepsilon)/n_\mathrm{co}\right]$ \ (photons)")

    fname = fdir + 'fig_gap_colorbar.pdf'
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    print("saved", fname)
