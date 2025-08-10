#-----------------------------------------------------------------------------#
#																			  #
# Created    : 2022															  #
# Author(s)  : Abdelrahman Hussein a.h.a.hussein@outlook.com				  #
#																			  #
# This file is part of the software `PHIMATS RVE`							  #
#																			  #
#-----------------------------------------------------------------------------#

import numpy as np

from PHIMATS_RVE import PhaseField
from numba_kernels import *

class Initializations:
    """
    A class to manage functions for initializing grains.
    """
   
#-------------------------------------------------------------------------------------------------#
        
    def Single(Phase: PhaseField, GID: int):
        """
        Initializes initial grain (back ground phase).

        Args:
            Phase (PhaseField): PhaseField object.
            GID (int): Grain id. 
        """
        
        if Phase.is3D:
            Phase.PhaseFields[GID, :, :, :] = 1
        else :
            Phase.PhaseFields[GID, :, :] = 1
                
#-------------------------------------------------------------------------------------------------#
        
    def Mirror(Phase: PhaseField, GID: int, axis: str = "x"):
        """
        Adds a second grain mirrored along `x`, `y` or `z` axis. 

        Args:
            Phase (PhaseField): PhaseField object.
            GID (int): Grain id.
            axis (str, optional): Axis. Options are `x`, `y` or `z`. Defaults to `x`.

        Raises:
            ValueError: If axis not equal to `x`, `y` or `z`. 
        """
        
        Nx, Ny = Phase.Nx, Phase.Ny
        if Phase.is3D:
            Nz = Phase.Nz

        if axis == "y":
            if Phase.is3D:
                Phase.PhaseFields[GID, int(Ny/2):, :, :] = 1
            else:
                Phase.PhaseFields[GID, int(Ny/2):, :] = 1

        elif axis == "x":
            if Phase.is3D:
                Phase.PhaseFields[GID, :, int(Nx/2):, :] = 1
            else:
                Phase.PhaseFields[GID, :, int(Nx/2):] = 1

        elif axis == "z":
            if not Phase.is3D:
                raise ValueError("Z-axis mirroring is only available for 3D fields.")
            Phase.PhaseFields[GID, :, :, int(Nz/2):] = 1

        else:
            raise ValueError("Invalid axis. Allowed values are  'x', 'y', or 'z'.")

        # Update phase 0
        Phase.PhaseFields[0, ...] = Phase.PhaseFields[0, ...] - Phase.PhaseFields[GID, ...]
            
#-------------------------------------------------------------------------------------------------#
        
    def Circle(Phase: PhaseField, GID: int, x0: float, y0: float, Radius: float = 3):
        """
        Initializes a circular phase.

        Args:
            Phase (PhaseField): PhaseField object.
            GID (int): Grain id.
            x0 (float): x coordinate.
            y0 (float): y coordinate.
            Radius (float, optional): Initial seed radius. Defaults to 3.
        """

        Nx = Phase.Nx 
        Ny = Phase.Ny
        
        dx = Phase.dx
        
        r_nuclei = Radius*dx
        
        www = Phase.www
        aaa = Phase.aaa

        x_nuclei = x0
        y_nuclei = y0
        
        for m in range(Ny):  # loop through Ny
            for l in range(Nx):  # loop through Nx

                # PBCs
                if l > Nx-1: 
                    l = l - Nx
                if l < 0:
                    l = l + Nx
                if m > Ny-1:
                    m = m - Ny
                if m < 0:
                    m = m + Ny

                r = np.sqrt((l*dx-x_nuclei*dx)**2 + \
                            (m*dx-y_nuclei*dx)**2)  \
                            - r_nuclei

                tmp = np.sqrt(2.*www)/aaa*r
                phi_tmp = 0.5*(1.-np.sin(tmp))

                if tmp >= np.pi/2.:
                    phi_tmp = 0.

                if phi_tmp > 0: 
                    Phase.PhaseFields[GID,l,m] = phi_tmp            
                    Phase.PhaseFields[0,l,m] = Phase.PhaseFields[0,l,m] - Phase.PhaseFields[GID,l,m]  

            # end l
            
        # end m
    
# -------------------------------------------------------------------------------------------------#

    def Seeds(Phase: PhaseField, nSeeds: int):
        """
        Initializes random seeds

        Args:
            Phase (PhaseField): PhaseField object.
            nSeeds (int): number of seeds.
        """
        
        print("Initilizing seeds ....")
                
        if Phase.is3D:
            Seeds_3D(Phase.PhaseFields, Phase.dx, Phase.www, Phase.aaa)
        else:
            Seeds_2D(Phase.PhaseFields, Phase.dx, Phase.www, Phase.aaa, nSeeds)
            
        Phase.UpdateGrain()
        
        print(nSeeds, "seeds initialized.\n")

#-------------------------------------------------------------------------------------------------#
