
"""
/***************************************************************************
 BufferByPercentage
                                 A (former) QGIS plugin
 Buffer polygon features so the buffered area is a specified percentage of
 the original area
                              -------------------
        begin                : 2013-10-12
        copyright            : (C) 2020 by Juernjakob Dugge
        
 edited by Leonardo Barleta, 2021-12-02
 convert scripts from QGIS Plugin to standalone Python functions
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os
import geopandas as gpd
import shapely
from shapely.geometry import Polygon
from pathlib import Path

def bufferByPercentage(geometry, target_factor, outer=True):
    buffer = geometry.buffer(findBufferLength(geometry, target_factor))
    
    return buffer if outer == False else buffer.difference(geometry)

def findBufferLength(geometry, target_factor, segments=10):
    """Find the buffer length that scales a geometry by a certain factor."""
    area_unscaled = geometry.area
    
    width = abs(geometry.bounds[0] - geometry.bounds[2])
    height = abs(geometry.bounds[1] - geometry.bounds[3])
    
    buffer_initial = 0.1 * (width + height)
    
    buffer_length = secant(calculateError, buffer_initial,
                           2 * buffer_initial, geometry, segments,
                           area_unscaled, target_factor)

    return buffer_length


def calculateError(buffer_length, geometry, segments, area_unscaled,
                   target_factor):
    """Calculate the difference between the current and the target factor."""
    geometry_scaled = geometry.buffer(buffer_length)
    area_scaled = geometry_scaled.area
    
    if area_scaled == 0:
        raise ValueError('Buffer length leads to zero-area polygon')

    return area_scaled / area_unscaled - target_factor


# Secant method for iteratively finding the root of a function
# Taken from
# http://www.physics.rutgers.edu/~masud/computing/WPark_recipes_in_python.html
def secant(func, oldx, x, *args, **kwargs):
    """Find the root of a function"""
    tolerance = kwargs.pop('tolerance', 1e-6)
    max_steps = kwargs.pop('max_steps', 100)

    steps = 0
    dx = 0
    oldf, f = func(oldx, *args), func(x, *args)

    if (abs(f) > abs(oldf)):  # Determine the initial search direction
        oldx, x = x, oldx
        oldf, f = f, oldf

    while (f - oldf) != 0 and steps < max_steps:
        dx = f * (x - oldx) / float(f - oldf)

        if abs(dx) < tolerance * (1 + abs(x)):  # Converged
            return x - dx

        oldx, x = x, x - dx
        
        try:
            oldf, f = f, func(x, *args)
        except ValueError:
            # The current step caused an invalid result. Halve the step size
            x = oldx  # Undo current step
            f = oldf
            dx *= 0.5  # Halve the step size
            oldx, x = x, x - dx
            oldf, f = f, func(x, *args)

        steps += 1

    # Did not converge
    return x - dx
