#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
utils.py

Routines to read and process data from sputniPIC simulations.
"""

# Built-in imports
import os

# Third-party imports
import numpy as np
from vtk import vtkStructuredPointsReader
from vtk.util.numpy_support import vtk_to_numpy

__author__ = "Louis Richard"
__email__ = "louisr@irf.se"
__copyright__ = "Copyright 2026"
__license__ = "MIT"

__all__ = [
    "read_input",
    "read_scalar",
    "read_field",
    "calc_gradx",
    "calc_grady",
    "calc_gradz",
    "calc_curl",
]


def read_input(path):
    r"""Reads the input file of a sputniPIC simulation and returns a 
    dictionary with the parameters.

     Parameters
     ----------
     path : str
         Path to the input file of the sputniPIC simulation.

    Returns
    -------
    meta : dict
        Dictionary with the parameters of the sputniPIC simulation.
    """

    with open(path, "r") as f:
        for i in range(4):
            f.readline()

        # Species properties
        n_species = int(f.readline().split(" = ")[1])

        nop, qom, rho0 = [np.zeros(n_species) for _ in range(3)]

        for i in range(n_species):
            line = f.readline()
            nop[i] = int(line.split("  ")[0].split(" = ")[1].split("\t")[0])
            qom[i] = int(line.split("  ")[1].split(" = ")[1])

        f.readline()
        # Box properties
        box_size = [float(f.readline().split("= ")[1]) for _ in range(3)]
        box_grid = [int(f.readline().split("= ")[1]) for _ in range(3)]

        f.readline()
        # Time step and number of cycle
        dt = float(f.readline().split("= ")[1])
        nt = int(f.readline().split("= ")[1])

        f.readline()
        for i in range(n_species):
            rho0[i] = float(f.readline().split("= ")[1])

        # Current sheet thickness
        h0 = float(f.readline().split("= ")[1])

        # Initial magnetic field
        b0_xyz = [float(f.readline().split("= ")[1]) for _ in range(3)]

        f.readline()
        smooth = int(f.readline().split("= ")[1])

        gmres_err_tol, cg_err_tol, mover_err_tol = [
            float(f.readline().split("= ")[1]) for _ in range(3)
        ]

    meta = {
        "n_species": n_species,
        "nop": nop,
        "qom": qom,
        "rho0": rho0,
        "box_size": box_size,
        "box_grid": box_grid,
        "dt": dt,
        "nt": nt,
        "h0": h0,
        "b0": b0_xyz,
        "smooth": smooth,
        "gmres_err_tol": gmres_err_tol,
        "cg_err_tol": cg_err_tol,
        "mover_err_tol": mover_err_tol,
    }

    return meta


def read_scalar(name, time_step, path):
    r"""Reads a scalar variable from a sputniPIC simulation and returns it as a 
    numpy array.

     Parameters
     ----------
     name : str
         Name of the variable to read. It should be the same as the name of the 
         variable in the sputniPIC simulation.
    time_step : int
        Time step to read. It should be the same as the time step of the variable 
        in the sputniPIC simulation.
    path : str
        Path to the directory containing the vtk files of the sputniPIC simulation.

    Returns
    -------
    scalar : numpy.ndarray
        Scalar variable read from the sputniPIC simulation.
    """

    reader = vtkStructuredPointsReader()
    reader.SetFileName(os.path.join(path, "{}_{:d}.vtk".format(name, int(time_step))))
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    dim = data.GetDimensions()
    scalar = vtk_to_numpy(data.GetPointData().GetArray(name.replace("_", "")))
    scalar = np.reshape(scalar, (dim[1], dim[0], dim[2]))

    return scalar


def read_field(name, time_step, path):
    r"""Reads a vector variable from a sputniPIC simulation and returns it as a numpy 
    array.

    Parameters
    ----------
    name : str
        Name of the variable to read. It should be the same as the name of the 
        variable in the sputniPIC simulation.
    time_step : int
        Time step to read. It should be the same as the time step of the variable 
        in the sputniPIC simulation.
    path : str
        Path to the directory containing the vtk files of the sputniPIC simulation.

    Returns
    -------
    field : numpy.ndarray
        Vector variable read from the sputniPIC simulation.

    """
    reader = vtkStructuredPointsReader()
    reader.SetFileName(os.path.join(path, "{}_{:d}.vtk".format(name, int(time_step))))
    reader.ReadAllVectorsOn()
    reader.Update()
    data = reader.GetOutput()
    dim = data.GetDimensions()
    field = vtk_to_numpy(data.GetPointData().GetArray(name))
    field = np.reshape(field, (dim[1], dim[0], dim[2], 3))

    return field


def calc_gradx(inp=None, x=None, periodic=True, der_ord=1):
    r"""Calculates x gradient component with(out) periodical boundary conditions.

    Parameters
    ----------
    inp : xarray.DataArray
        Input variable of the procedure

    x : xarray.DataArray
        x axis coordinates.

    periodic : bool
        Periodic boundary conditions? Default is True.

    der_ord : int
        Order of derivation 
        (0: no derivative, 1: first derivative, 2: second derivative ...).
        Default is 1.

    Returns
    -------
    dx_f : numpy.ndarray
        Gradient of the input in the x direction

    """

    # determine the coefficients
    ny, nx, nz = inp.shape

    dx = np.median(np.diff(x.data))

    if nx == 1:
        return np.zeros([nx, ny, nz])

    else:
        oo_6 = dx ** (-der_ord) * np.array([1.0, 0.0, -49.0 / 18.0])
        aa_6 = dx ** (-der_ord) * np.array([0.0, 9.0 / 12.0, 3.0 / 2.0])
        bb_6 = dx ** (-der_ord) * np.array([0.0, -3.0 / 20.0, -3.0 / 20.0])
        cc_6 = dx ** (-der_ord) * np.array([0.0, 1.0 / 60.0, 1.0 / 90.0])

        # create the new vector, fill it
        if periodic:
            ff = np.tile(inp.data, (1, 2, 1))
            dx_f = oo_6[der_ord] * ff[:, 0:nx, :]
            dx_f += aa_6[der_ord] * (
                ff[:, 1 : 1 + nx, :] - ff[:, nx - 1 : 2 * nx - 1, :]
            )
            dx_f += bb_6[der_ord] * (
                ff[:, 2 : 2 + nx, :] - ff[:, nx - 2 : 2 * nx - 2, :]
            )
            dx_f += cc_6[der_ord] * (
                ff[:, 3 : 3 + nx, :] - ff[:, nx - 3 : 2 * nx - 3, :]
            )
        else:
            f = inp.data
            dx_f = np.zeros([ny, nx, nz])
            dx_f[:, 3 : nx - 3, :] = oo_6[der_ord] * f[:, 3 : nx - 3, :]
            dx_f[:, 3 : nx - 3, :] += aa_6[der_ord] * (
                f[:, 4 : nx - 2, :] - f[:, 2 : nx - 4, :]
            )
            dx_f[:, 3 : nx - 3, :] += bb_6[der_ord] * (
                f[:, 5 : nx - 1, :] - f[:, 1 : nx - 5, :]
            )
            dx_f[:, 3 : nx - 3, :] += cc_6[der_ord] * (
                f[:, 6 : nx - 0, :] - f[:, 0 : nx - 6, :]
            )

        # print('done with the new array!')
        return dx_f


def calc_grady(inp=None, y=None, periodic=True, der_ord=1):
    r"""Calculates y gradient component with(out) periodical boundary conditions.

    Parameters
    ----------
    inp : xarray.DataArray
        Input variable of the procedure

    y : xarray.DataArray
        y axis coordinates.

    periodic : bool
        Periodic boundary conditions? Default is True.

    der_ord : int
        Order of derivation 
        (0: no derivative, 1: first derivative, 2: second derivative ...).
        Default is 1.

    Returns
    -------
    dy_f : numpy.ndarray
        Gradient of the input in the y direction

    """

    # determine the coefficients
    ny, nx, nz = inp.shape

    dy = np.median(np.diff(y.data))

    if ny == 1:
        return np.zeros([ny, nx, nz])

    else:
        oo_6 = dy ** (-der_ord) * np.array([1.0, 0.0, -49.0 / 18.0])
        aa_6 = dy ** (-der_ord) * np.array([0.0, 9.0 / 12.0, 3.0 / 2.0])
        bb_6 = dy ** (-der_ord) * np.array([0.0, -3.0 / 20.0, -3.0 / 20.0])
        cc_6 = dy ** (-der_ord) * np.array([0.0, 1.0 / 60.0, 1.0 / 90.0])

        # create the new vector, fill it
        if periodic:
            ff = np.tile(inp.data, (2, 1, 1))
            dy_f = oo_6[der_ord] * ff[0:ny, :, :]
            dy_f += aa_6[der_ord] * (
                ff[1 : 1 + ny, :, :] - ff[ny - 1 : 2 * ny - 1, :, :]
            )
            dy_f += bb_6[der_ord] * (
                ff[2 : 2 + ny, :, :] - ff[ny - 2 : 2 * ny - 2, :, :]
            )
            dy_f += cc_6[der_ord] * (
                ff[3 : 3 + ny, :, :] - ff[ny - 3 : 2 * ny - 3, :, :]
            )
        else:
            f = inp.data
            dy_f = np.zeros([ny, nx, nz])
            dy_f[3 : ny - 3, :, :] = oo_6[der_ord] * f[3 : ny - 3, :, :]
            dy_f[3 : ny - 3, :, :] += aa_6[der_ord] * (
                f[4 : ny - 2, :, :] - f[2 : ny - 4, :, :]
            )
            dy_f[3 : ny - 3, :, :] += bb_6[der_ord] * (
                f[5 : ny - 1, :, :] - f[1 : ny - 5, :, :]
            )
            dy_f[3 : ny - 3, :, :] += cc_6[der_ord] * (
                f[6 : ny - 0, :, :] - f[0 : ny - 6, :, :]
            )

        # print('done with the new array!')
        return dy_f


def calc_gradz(inp=None, z=None, periodic=True, der_ord=1):
    r"""Calculates z gradient component with(out) periodical boundary conditions.

    Parameters
    ----------
    inp : xarray.DataArray
        Input variable of the procedure

    z : xarray.DataArray
        z axis coordinates.

    periodic : bool
        Periodic boundary conditions? Default is True.

    der_ord : int
        Order of derivation 
        (0: no derivative, 1: first derivative, 2: second derivative ...).
        Default is 1.

    Returns
    -------
    dz_f : numpy.ndarray
        Gradient of the input in the y direction

    """

    # determine the coefficients
    ny, nx, nz = inp.shape

    dz = np.median(np.diff(z.data))

    if nz == 1:
        return np.zeros([ny, nx, nz])

    else:
        oo_6 = dz ** (-der_ord) * np.array([1.0, 0.0, -49.0 / 18.0])
        aa_6 = dz ** (-der_ord) * np.array([0.0, 9.0 / 12.0, 3.0 / 2.0])
        bb_6 = dz ** (-der_ord) * np.array([0.0, -3.0 / 20.0, -3.0 / 20.0])
        cc_6 = dz ** (-der_ord) * np.array([0.0, 1.0 / 60.0, 1.0 / 90.0])

        # create the new vector, fill it
        if periodic:
            ff = np.tile(inp.data, (1, 1, 2))
            dz_f = oo_6[der_ord] * ff[:, :, 0:nz]
            dz_f += aa_6[der_ord] * (
                ff[:, :, 1 : 1 + nz] - ff[:, :, nz - 1 : 2 * nz - 1]
            )
            dz_f += bb_6[der_ord] * (
                ff[:, :, 2 : 2 + nz] - ff[:, :, nz - 2 : 2 * nz - 2]
            )
            dz_f += cc_6[der_ord] * (
                ff[:, :, 3 : 3 + nz] - ff[:, :, nz - 3 : 2 * nz - 3]
            )
        else:
            f = inp.data
            dz_f = np.zeros([ny, nx, nz])
            dz_f[:, :, 3 : nz - 3] = oo_6[der_ord] * f[:, :, 3 : nz - 3]
            dz_f[:, :, 3 : nz - 3] += aa_6[der_ord] * (
                f[:, :, 4 : nz - 2] - f[:, :, 2 : nz - 4]
            )
            dz_f[:, :, 3 : nz - 3] += bb_6[der_ord] * (
                f[:, :, 5 : nz - 1] - f[:, :, 1 : nz - 5]
            )
            dz_f[:, :, 3 : nz - 3] += cc_6[der_ord] * (
                f[:, :, 6 : nz - 0] - f[:, :, 0 : nz - 6]
            )

        # print('done with the new array!')
        return dz_f


def calc_curl(inp, x=None, y=None, z=None, perx=True, pery=False, perz=True):
    r"""Calculates curl of the input field with(out) periodical boundary conditions.

    Parameters
    ----------
    inp : xarray.DataArray
        Input variable of the procedure

    x : xarray.DataArray
        x axis coordinates.

    y : xarray.DataArray
        y axis coordinates.

    z : xarray.DataArray
        x axis coordinates.

    perx : bool
        Periodic boundary conditions in the x direction ? Delault is True.

    pery : bool
        Periodic boundary conditions in the y direction ? Delault is True.

    perz : bool
        Periodic boundary conditions in the z direction ? Delault is True.

    Returns
    -------
    dy_f : numpy.ndarray
        Gradient of the input in the y direction

    """
    # create a raw vector field
    ny, nx, nz = inp.shape[:-1]
    curl_x = np.zeros([ny, nx, nz])
    curl_y = np.zeros([ny, nx, nz])
    curl_z = np.zeros([ny, nx, nz])

    if nx != 1:
        curl_y -= calc_gradx(inp[..., 2], x, periodic=perx)
        curl_z += calc_gradx(inp[..., 1], x, periodic=perx)
    if ny != 1:
        curl_z -= calc_grady(inp[..., 0], y, periodic=pery)
        curl_x += calc_grady(inp[..., 2], y, periodic=pery)
    if nz != 1:
        curl_x -= calc_gradz(inp[..., 1], z, periodic=perz)
        curl_y += calc_gradz(inp[..., 0], z, periodic=perz)

    out = np.zeros(inp.shape)
    out[..., 0] = curl_x
    out[..., 1] = curl_y
    out[..., 2] = curl_z

    return out
