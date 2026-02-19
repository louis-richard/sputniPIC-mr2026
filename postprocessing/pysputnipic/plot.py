#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot.py

Routines to plot data from sputniPIC simulations.
"""

import matplotlib.pyplot as plt
import numpy as np

__author__ = "Louis Richard"
__email__ = "louisr@irf.se"
__copyright__ = "Copyright 2026"
__license__ = "MIT"

__all__ = ["plot_all_2d", "plot_2d", "plot_spacecraft"]


def plot_all_2d(inp=None, time=0):
    # Get limits for all variables
    # Symmetric limits for B, E, ve, and vi
    bmin = [np.min(inp.b_xyz.isel(z=0, comp=i).data) for i in range(3)]
    bmax = [np.max(inp.b_xyz.isel(z=0, comp=i).data) for i in range(3)]

    bmin[0], bmax[0] = [np.max([np.abs(bmin[0]), bmax[0]]) for _ in range(2)]
    bmin[1], bmax[1] = [np.max([np.abs(bmin[1]), bmax[1]]) for _ in range(2)]
    bmin[2], bmax[2] = [np.max([np.abs(bmin[2]), bmax[2]]) for _ in range(2)]
    bmin = -np.array(bmin)

    emin = [np.min(inp.e_xyz.isel(z=0, comp=i).data) for i in range(3)]
    emax = [np.max(inp.e_xyz.isel(z=0, comp=i).data) for i in range(3)]

    emin[0], emax[0] = [np.max([np.abs(emin[0]), emax[0]]) for _ in range(2)]
    emin[1], emax[1] = [np.max([np.abs(emin[1]), emax[1]]) for _ in range(2)]
    emin[2], emax[2] = [np.max([np.abs(emin[2]), emax[2]]) for _ in range(2)]
    emin = -np.array(emin)

    vemin = [np.min(inp.ve_xyz.isel(z=0, comp=i).data) for i in range(3)]
    vemax = [np.max(inp.ve_xyz.isel(z=0, comp=i).data) for i in range(3)]

    vemin[0], vemax[0] = [np.max([np.abs(vemin[0]), vemax[0]]) for _ in range(2)]
    vemin[1], vemax[1] = [np.max([np.abs(vemin[1]), vemax[1]]) for _ in range(2)]
    vemin[2], vemax[2] = [np.max([np.abs(vemin[2]), vemax[2]]) for _ in range(2)]
    vemin = -np.array(vemin)

    vimin = [np.min(inp.vi_xyz.isel(z=0, comp=i).data) for i in range(3)]
    vimax = [np.max(inp.vi_xyz.isel(z=0, comp=i).data) for i in range(3)]

    vimin[0], vimax[0] = [np.max([np.abs(vimin[0]), vimax[0]]) for _ in range(2)]
    vimin[1], vimax[1] = [np.max([np.abs(vimin[1]), vimax[1]]) for _ in range(2)]
    vimin[2], vimax[2] = [np.max([np.abs(vimin[2]), vimax[2]]) for _ in range(2)]
    vimin = -np.array(vimin)

    rhomin = [np.min(inp.rho.isel(z=0, specie=i).data) for i in range(3)]
    rhomax = [np.max(inp.rho.isel(z=0, specie=i).data) for i in range(3)]

    rhomin[2], rhomax[2] = [np.max([np.abs(rhomin[2]), rhomax[2]]) for _ in range(2)]
    rhomin = np.array([rhomin[0], rhomin[1], -rhomin[2]])

    # Plot
    fig, axs = plt.subplots(5, 3, figsize=(16, 14))
    fig.subplots_adjust(
        bottom=0.05, top=0.95, left=0.05, right=0.95, hspace=0.4, wspace=0.3
    )

    for i in range(3):
        # b
        kwargs_b = {
            "vmin": bmin[i],
            "vmax": bmax[i],
            "interpolation": "bicubic",
            "cmap": "RdBu",
        }
        inp.b_xyz.isel(time=time, z=0, comp=i).plot.imshow(
            ax=axs[0, i], x="x", y="y", **kwargs_b
        )

        # e
        kwargs_e = {
            "vmin": emin[i],
            "vmax": emax[i],
            "interpolation": "bicubic",
            "cmap": "RdBu",
        }
        inp.e_xyz.isel(time=time, z=0, comp=i).plot.imshow(
            ax=axs[1, i], x="x", y="y", **kwargs_e
        )

        # ve
        kwargs_ve = {
            "vmin": vemin[i],
            "vmax": vemax[i],
            "interpolation": "bicubic",
            "cmap": "RdBu",
        }
        inp.ve_xyz.isel(time=time, z=0, comp=i).plot.imshow(
            ax=axs[2, i], x="x", y="y", **kwargs_ve
        )

        # vi
        kwargs_vi = {
            "vmin": vimin[i],
            "vmax": vimax[i],
            "interpolation": "bicubic",
            "cmap": "RdBu",
        }
        inp.vi_xyz.isel(time=time, z=0, comp=i).plot.imshow(
            ax=axs[3, i], x="x", y="y", **kwargs_vi
        )

    # rho
    kwargs_rho_i = {
        "vmin": rhomin[0],
        "vmax": rhomax[0],
        "interpolation": "bicubic",
        "cmap": "Blues",
    }
    inp.rho.isel(time=time, z=0, specie=0).plot.imshow(
        ax=axs[4, 0], x="x", y="y", **kwargs_rho_i
    )

    kwargs_rho_e = {
        "vmin": rhomin[1],
        "vmax": rhomax[1],
        "interpolation": "bicubic",
        "cmap": "Reds_r",
    }
    inp.rho.isel(time=time, z=0, specie=1).plot.imshow(
        ax=axs[4, 1], x="x", y="y", **kwargs_rho_e
    )

    kwargs_rho_n = {
        "vmin": rhomin[2],
        "vmax": rhomax[2],
        "interpolation": "bicubic",
        "cmap": "RdBu",
    }
    inp.rho.isel(time=time, z=0, specie=2).plot.imshow(
        ax=axs[4, 2], x="x", y="y", **kwargs_rho_n
    )

    return fig, axs


def plot_2d(inp=None, tar_var="", time=0, scx=None, show_sc=True):
    x = inp["b_xyz"].x.data
    y = inp["b_xyz"].y.data
    b_x = np.squeeze(inp["b_xyz"].data[time, :, :, :, 0])
    b_y = np.squeeze(inp["b_xyz"].data[time, :, :, :, 1])

    if tar_var == "rho":
        g = (
            inp[tar_var]
            .isel(time=time, z=0)
            .plot.imshow("x", "y", col="specie", col_wrap=3, interpolation="bicubic")
        )

    else:
        g = (
            inp[tar_var]
            .isel(time=time, z=0)
            .plot.imshow("x", "y", col="comp", col_wrap=3, interpolation="bicubic")
        )

    for ax in list(g.axes[0, :]):
        ax.streamplot(
            x,
            y,
            b_x,
            b_y,
            color="k",
            linewidth=1,
            density=1,
            arrowstyle="->",
            arrowsize=0.8,
        )

    if scx is not None and isinstance(scx, int):
        for ax in list(g.axes[0, :]):
            ax.axvline(inp.x.data[scx], linestyle="--", color="k")

        f, axs = plot_spacecraft(inp, scx, time)

        if not show_sc:
            plt.close(f)

    return


def plot_spacecraft(inp=None, time=0, x=64):
    f, axs = plt.subplots(4, sharex="all", figsize=(6, 9))
    f.subplots_adjust(left=0.15, right=0.85, bottom=0.05, top=0.95, hspace=0)
    inp.b_xyz.isel(x=64, z=0, comp=[0, 1, 2], time=time).plot.line(ax=axs[0], x="y")
    axs[0].set_xlabel("")
    inp.e_xyz.isel(x=64, z=0, comp=[0, 1, 2], time=time).plot.line(ax=axs[1], x="y")
    axs[1].set_xlabel("")
    axs[1].set_title("")
    inp.j_xyz.isel(x=64, z=0, comp=[0, 1, 2], time=time).plot.line(ax=axs[2], x="y")
    axs[2].set_xlabel("")
    axs[2].set_title("")
    inp.rho.isel(x=64, z=0, specie=[0, 1, 2], time=time).plot.line(ax=axs[3], x="y")
    axs[3].set_title("")

    f.align_ylabels(axs)

    return f, axs
