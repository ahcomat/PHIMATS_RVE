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
        initializes initial grain (back ground phase).

        Args:
            Phase (PhaseField): PhaseField object.
            GID (int): Grain id. 
        """
        
        if Phase.is3D:
            Phase.PhaseFields[GID, :, :, :] = 1
        else :
            Phase.PhaseFields[GID, :, :] = 1
                
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
