#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gem_2d.py

Script to read and plot 2D GEM data from sputniPIC simulations.
"""

import argparse
import glob

# Built-in imports
import os

import matplotlib.pyplot as plt

# Third-party imports
import numpy as np
import xarray as xr
from tqdm import tqdm
from scipy.ndimage import gaussian_filter

# Local imports
from pysputnipic import utils
from pysputnipic.plot import plot_all_2d

__author__ = "Louis Richard"
__email__ = "louisr@irf.se"
__copyright__ = "Copyright 2026"
__license__ = "MIT"


def main(args):
    if os.path.isfile(args.path):
        gem_2d = xr.load_dataset(args.path)
    else:
        path, _ = os.path.split(args.path)

        # Get times
        times = [
            int(os.path.split(f)[1].strip(".vtk").split("_")[1])
            for f in glob.glob("{}/B_*".format(path))
        ]
        times = np.sort(times)
        print(times)

        if args.time is not None:
            if args.time in times:
                times = np.array([args.time])
            else:
                raise ValueError("Time step {} not found in data".format(args.time))

        # Number of time steps to read and plot
        nt = len(times)

        # Name species
        species = ["ion", "electron", "net"]

        # Get box size
        # Read meta
        meta = utils.read_input(os.path.join(path, "sputniPICparameters.txt"))

        # Time step
        dt = meta["dt"]

        # Boz size
        lx, ly, lz = meta["box_size"]
        nx, ny, nz = meta["box_grid"]

        x, y, z = [
            np.linspace(0, l, n) for l, n in zip(meta["box_size"], meta["box_grid"])
        ]

        # Initialize arrays
        b_xyz, e_xyz = [np.zeros((nt, ny, nx, nz, 3)) for _ in range(2)]
        ve_xyz, vi_xyz = [np.zeros((nt, ny, nx, nz, 3)) for _ in range(2)]

        rho_i, rho_e, rho_n = [np.zeros((nt, ny, nx, nz)) for _ in range(3)]

        for it, t in enumerate(tqdm(times)):
            # Read electromagnetic fields
            b_xyz[it, ...] = utils.read_field("B", t, path)
            e_xyz[it, ...] = utils.read_field("E", t, path)

            # Read electron and ion velocities
            ve_xyz[it, ...] = utils.read_field("Ve", t, path)
            vi_xyz[it, ...] = utils.read_field("Vi", t, path)

            # Read densities
            rho_i[it, ...] = utils.read_scalar("rho{}".format("i"), t, path)
            rho_e[it, ...] = utils.read_scalar("rho{}".format("e"), t, path)
            rho_n[it, ...] = utils.read_scalar("rho{}".format("_net"), t, path)

        # Filter 
        sigma = 2
        b_xyz = gaussian_filter(b_xyz, sigma=sigma, axes=(1, 2))
        e_xyz = gaussian_filter(e_xyz, sigma=sigma, axes=(1, 2))
        ve_xyz = gaussian_filter(ve_xyz, sigma=sigma, axes=(1, 2))
        vi_xyz = gaussian_filter(vi_xyz, sigma=sigma, axes=(1, 2))
        rho_i = gaussian_filter(rho_i, sigma=sigma, axes=(1, 2))
        rho_e = gaussian_filter(rho_e, sigma=sigma, axes=(1, 2))
        rho_n = gaussian_filter(rho_n, sigma=sigma, axes=(1, 2))

        # Convert to xarray
        b_xyz = xr.DataArray(
            b_xyz,
            coords=[times * dt, y, x, z, ["x", "y", "z"]],
            dims=["time", "y", "x", "z", "comp"],
        )
        e_xyz = xr.DataArray(
            e_xyz,
            coords=[times * dt, y, x, z, ["x", "y", "z"]],
            dims=["time", "y", "x", "z", "comp"],
        )

        ve_xyz = xr.DataArray(
            ve_xyz,
            coords=[times * dt, y, x, z, ["x", "y", "z"]],
            dims=["time", "y", "x", "z", "comp"],
        )
        vi_xyz = xr.DataArray(
            vi_xyz,
            coords=[times * dt, y, x, z, ["x", "y", "z"]],
            dims=["time", "y", "x", "z", "comp"],
        )

        rho = np.stack([rho_i, rho_e, rho_n])
        rho = xr.DataArray(
            rho,
            coords=[species, times * dt, y, x, z],
            dims=["specie", "time", "y", "x", "z"],
        )

        gem_2d = xr.Dataset(
            {
                "b_xyz": b_xyz,
                "e_xyz": e_xyz,
                "ve_xyz": ve_xyz,
                "vi_xyz": vi_xyz,
                "rho": rho,
            }
        )

        # gem_2d.to_netcdf(args.path)

    if args.movie:
        for it in tqdm(np.arange(len(gem_2d.time))):
            fig, _ = plot_all_2d(gem_2d, it)
            fig_path = os.path.join(os.path.split(path.strip("/"))[0], "figures")
            fig.savefig(os.path.join(fig_path, "image_{:03d}.png".format(times[it])))
            fig.close()
    else:
        it = 0
        fig, _ = plot_all_2d(gem_2d, it)
        fig_path = os.path.join(os.path.split(path.strip("/"))[0], "figures")
        print(fig_path)
        fig.savefig(os.path.join(fig_path, "image_{:03d}.png".format(times[it])))
        plt.show()

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="data path", action="store")
    parser.add_argument("-t", "--time", help="time step", action="store", type=int)
    parser.add_argument("-m", "--movie", help="time step", action="store_true")

    main(parser.parse_args())
