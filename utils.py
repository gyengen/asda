#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 2019

Name: utils.py

Purpose: Util functions for ASDA

@author: Jaijia Liu at University of Sheffield
"""
__author__ = 'Jiajia Liu'
__copyright__ = 'Copyright 2019, The Solar Physics and Space Plasma ' + \
                'Research Center (SP2RC)'
__license__ = 'GPLv3'
__date__ = '2019/11/27'
__maintainor__ = 'Jiajia Liu'
__email__ = 'jj.liu@sheffield.ac.uk'


import numpy as np
from skimage import measure
from scipy.interpolate import interp2d


def reform2d(array, factor=1):
    """
    Reform a 2d array by a given factor

    Parameters
    ----------
    array : `numpy.ndarray`
        2d array to be reformed

    factor : `int`
        The array is going to be magnified by the factor. Default is 1.

    Returns
    -------
        `numpy.ndarray`
        reformed array

    """
    if not isinstance(factor, int):
        raise ValueError("Parameter 'factor' must be an integer!")

    if len(np.shape(array)) != 2:
        raise ValueError("Input array must be 2d!")

    if factor > 1:
        congridx = interp2d(np.arange(0, array.shape[0]),
                            np.arange(0, array.shape[1]), array.T)
        array = congridx(np.arange(0, array.shape[0], 1/factor),
                         np.arange(0, array.shape[1], 1/factor)).T

    return array


def points_in_poly(poly):
    """
    Return polygon as grid of points inside polygon. Only works for polygons
    defined with points which are all integers

    Parameters
    ----------
    poly : `list` or `numpy.ndarray`
        n x 2 list, defines all points at the edge of a polygon

    Returns
    -------
        `list`
        n x 2 array, all points within the polygon

    """
    if np.shape(poly)[1] != 2:
        raise ValueError("Polygon must be defined as a n x 2 array!")

    # convert to integers
    poly = np.array(poly, dtype=int).tolist()

    xs, ys = zip(*poly)
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    # New polygon with the staring point as [0, 0]
    newPoly = [(int(x - minx), int(y - miny)) for (x, y) in poly]
    mask = measure.grid_points_in_poly((round(maxx - minx) + 1,
                                        round(maxy - miny) + 1), newPoly)
    # all points in polygon
    points = [[x + minx, y + miny] for (x, y) in zip(*np.nonzero(mask))]

    # add edge points if missing
    for p in poly:
        if p not in points:
            points.append(p)

    return points
