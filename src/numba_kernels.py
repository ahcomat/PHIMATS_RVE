#-----------------------------------------------------------------------------#
#																			  #
# Created    : 2022															  #
# Author(s)  : Abdelrahman Hussein a.h.a.hussein@outlook.com				  #
#																			  #
# This file is part of the software `PHIMATS RVE`							  #
#																			  #
#-----------------------------------------------------------------------------#

import numpy as np
from numba import njit

#-------------------------------------------------------------------------------------------------#

# For PhaseField class

#-------------------------------------------------------------------------------------------------#

@njit()
def n_Calc_gPhi_2D(Field: np.ndarray) -> np.ndarray: 
    """
    Calculates the interface function gPhi for a 2D phase-field array.

    Args:
        Field (ndarray): Array of phase-field values.

    Returns:
        gPhi (ndarray): Interface array. 
    """
    NGrains, Nx, Ny = Field.shape

    gPhi = np.zeros((Nx, Ny))

    for i in range(NGrains):
        for j in range(NGrains):

            if i>j:
                gPhi += Field[i]*Field[j]

    return gPhi  

#-------------------------------------------------------------------------------------------------#

@njit()
def n_Calc_gPhi_3D(Field: np.ndarray) -> np.ndarray:
    """
    Computes the interface function gPhi for a 3D phase-field array.

    Args:
        Field (ndarray): Array of phase-field values.

    Returns:
        gPhi (ndarray): Diffuse interface array. 
    """
    
    NGrains, Nx, Ny, Nz = Field.shape

    gPhi = np.zeros((Nx, Ny, Nz))

    for i in range(NGrains):
        for j in range(NGrains):

            if i>j:
                gPhi += Field[i]*Field[j]

    return gPhi

#-------------------------------------------------------------------------------------------------#

@njit()
def n_UpdateGrain_2D(GId: np.ndarray, GNum: np.ndarray, Field: np.ndarray):
    """
    Updates grain IDs

    Args:
        GId (np.ndarray): _description_
        GNum (np.ndarray): _description_
        Field (np.ndarray): _description_
    """
                
    number_of_grain, Nx, Ny = Field.shape

    for m in range(Ny):
        for l in range(Nx):

            l_p = l + 1
            l_m = l - 1
            m_p = m + 1
            m_m = m - 1
            if l_p > Nx-1:
                l_p = l_p - Nx
            if l_m < 0:
                l_m = l_m + Nx
            if m_p > Ny-1:
                m_p = m_p - Ny
            if m_m < 0:
                m_m = m_m + Ny

            n = 0

            for i in range(number_of_grain):

                if Field[i,l,m] > 0.0 or (Field[i,l,m] == 0.0 and Field[i,l_p,m] > 0.0 or Field[i,l_m,m] > 0.0 or Field[i,l,m_p] > 0.0 or Field[i,l,m_m] > 0.0):

                    n += 1
                    GId[n-1,l,m] = i
                    GNum[l,m] = n

#-------------------------------------------------------------------------------------------------#

@njit()
def n_UpdateGrain_3D(GId, GNum, Field):
    """
    Updates grain IDs

    Args:
        GId (np.ndarray): _description_
        GNum (np.ndarray): _description_
        Field (np.ndarray): _description_
    """
    
    number_of_grain, Nx, Ny, Nz = Field.shape
        
    for k in range(Nz):
        for m in range(Ny):
            for l in range(Nx):

                l_p = l + 1
                l_m = l - 1
                m_p = m + 1
                m_m = m - 1
                k_p = k + 1
                k_m = k - 1
                
                if l_p > Nx-1:
                    l_p = l_p - Nx
                if l_m < 0:
                    l_m = l_m + Nx
                if m_p > Ny-1:
                    m_p = m_p - Ny
                if m_m < 0:
                    m_m = m_m + Ny
                if k_p > Nz-1:
                    k_p = k_p - Nz
                if k_m < 0:
                    k_m = k_m + Nz

                n = 0

                for i in range(number_of_grain):

                    if Field[i,l,m,k] > 0.0 or (Field[i,l,m,k] == 0.0 and Field[i,l_p,m,k] > 0.0 \
                        or Field[i,l_m,m,k] > 0.0 or Field[i,l,m_p,k] > 0.0 or Field[i,l,m_m,k] > 0.0 \
                            or Field[i,l,m,k_p] > 0.0 or Field[i,l,m,k_m] > 0.0):

                        n += 1
                        GId[n-1,l,m,k] = i

                GNum[l,m,k] = n

#-------------------------------------------------------------------------------------------------#

@njit()
def n_CalcTimeStep_2D(Field, dx, dt, PhiNew, GId, GNum, eij, wij, mij, aij):
    
    number_of_grain, Nx, Ny = Field.shape

    for m in range(Ny):
        for l in range(Nx):

            l_p1 = l + 1
            l_m1 = l - 1

            m_p1 = m + 1
            m_m1 = m - 1

            if l_p1 > Nx-1: 
                l_p1 = l_p1 - Nx

            if m_p1 > Ny-1: 
                m_p1 = m_p1 - Ny

            if l_m1 < 0:
                l_m1 = l_m1 + Nx

            if m_m1 < 0:
                m_m1 = m_m1 + Nx

            for n1 in range(GNum[l,m]):

                i = GId[n1, l, m]
                dpi = 0.0

                for n2 in range(GNum[l,m]):

                    j = GId[n2, l, m]
                    ppp = 0.0
                    Phii_Phij = Field[i,l,m]*Field[j,l,m]

                    for n3 in range(GNum[l,m]):

                        k = GId[n3, l, m]
                                    
                        laplacian = + Field[k, l_m1, m] \
                                    + Field[k, l, m_m1] \
                                    - 4*Field[k, l, m]  \
                                    + Field[k, l, m_p1] \
                                    + Field[k, l_p1, m] 

                        
                        ppp +=   (wij[i,k]-wij[j,k])*Field[k,l,m] \
                               + 0.5*(aij[i,k]**2 - aij[j,k]**2)*laplacian/dx**2

                    dpi = dpi - 2.0 * mij[i,j]/GNum[l,m] * (ppp - 8./np.pi*np.sqrt(Phii_Phij)*eij[i,j])

                    # end n3

                PhiNew[i,l,m] = Field[i,l,m] + dpi*dt

                # end n2
                
            # end n1

    # ----------<<< Finalize >>>----------    

    PhiNew = np.where(PhiNew <= 0.0,0.0, PhiNew)
    PhiNew = np.where(PhiNew >= 1.0,1.0, PhiNew)

    for m in range(Ny):
        for l in range(Nx):
            a = np.sum(PhiNew[:,l,m])
            Field[:,l,m] = PhiNew[:,l,m]/a

#-------------------------------------------------------------------------------------------------#

@njit()
def n_CalcTimeStep_3D(Field, dx, dt, PhiNew, GId, GNum, eij, wij, mij, aij):
    
    number_of_grain, Nx, Ny, Nz = Field.shape

    for n in range(Nz):
        for m in range(Ny):
            for l in range(Nx):

                l_p1 = l + 1
                l_m1 = l - 1

                m_p1 = m + 1
                m_m1 = m - 1
                
                n_p1 = n + 1
                n_m1 = n - 1


                if l_p1 > Nx-1: 
                    l_p1 = l_p1 - Nx

                if m_p1 > Ny-1: 
                    m_p1 = m_p1 - Ny
                    
                if n_p1 > Nz-1: 
                    n_p1 = n_p1 - Nz

                if l_m1 < 0:
                    l_m1 = l_m1 + Nx

                if m_m1 < 0:
                    m_m1 = m_m1 + Ny

                if n_m1 < 0:
                    n_m1 = n_m1 + Nz
                    
                for n1 in range(GNum[l,m,n]):

                    i = GId[n1, l, m, n]
                    dpi = 0.0

                    for n2 in range(GNum[l,m,n]):

                        j = GId[n2, l, m, n]
                        ppp = 0.0
                        Phii_Phij = Field[i,l,m,n]*Field[j,l,m,n]

                        for n3 in range(GNum[l,m,n]):

                            k = GId[n3, l, m, n]
                                        
                            laplacian = + Field[k, l_m1, m, n] \
                                        + Field[k, l, m_m1, n] \
                                        + Field[k, l, m, n_m1] \
                                        - 6*Field[k, l, m, n]  \
                                        + Field[k, l, m, n_p1] \
                                        + Field[k, l, m_p1, n] \
                                        + Field[k, l_p1, m, n] 
                            
                            ppp +=   (wij[i,k]-wij[j,k])*Field[k,l,m,n] \
                                + 0.5*(aij[i,k]**2 - aij[j,k]**2)*laplacian/dx**2

                        dpi = dpi - 2.0 * mij[i,j]/GNum[l,m,n] * (ppp - 8./np.pi*np.sqrt(Phii_Phij)*eij[i,j])

                        # end n3

                    PhiNew[i,l,m,n] = Field[i,l,m,n] + dpi*dt
                    
                    # end n2
                    
                # end n1

    # ----------<<< Finalize >>>----------    

    PhiNew = np.where(PhiNew <= 0.0,0.0, PhiNew)
    PhiNew = np.where(PhiNew >= 1.0,1.0, PhiNew)
    for n in range(Nz):
        for m in range(Ny):
            for l in range(Nx):
                a = np.sum(PhiNew[:,l,m,n])
                Field[:,l,m,n] = PhiNew[:,l,m,n]/a

#-------------------------------------------------------------------------------------------------#

@njit()
def n_getPFID_3D(Field: np.ndarray, OPid: np.ndarray, mf: np.ndarray, nf: np.ndarray):
    """
    Evaluates the grain ID for the order parameters.

    Args:
        Field (np.ndarray): Order parameters array.
        OPid (np.ndarray): Array for order parameters ids.
        mf (np.ndarray): _description_
        nf (np.ndarray): _description_
    """

    NumGrains, Nx, Ny, Nz = Field.shape

    for k in range(0, Nz):
        for m in range(0, Ny):
            for l in range(0, Nx):

                num_grain = 0
                max_phi = 0.

                for n in range(nf[l,m,k]):

                    n1 = mf[n,l,m,k]

                    if Field[n1,l,m,k] > max_phi:

                        max_phi = Field[n1,l,m,k]
                        num_grain = n1
                        OPid[l,m,k] = n1
                        
#-------------------------------------------------------------------------------------------------#
                        
@njit()
def n_getPFID_2D(Field: np.ndarray, OPid: np.ndarray, mf: np.ndarray, nf: np.ndarray):
    """
    Evaluates the grain ID for the order parameters.

    Args:
        Field (np.ndarray): Order parameters array.
        OPid (np.ndarray): Array for order parameters ids.
        mf (np.ndarray): _description_
        nf (np.ndarray): _description_
    """
        
    NumGrains, Nx, Ny = Field.shape

    for m in range(0, Ny):
        for l in range(0, Nx):

            max_phi = 0.

            for n in range(nf[l,m]):

                n1 = mf[n,l,m]

                if Field[n1,l,m] > max_phi:

                    max_phi = Field[n1,l,m]
                    OPid[l,m] = n1
        
#-------------------------------------------------------------------------------------------------#

# For Initialization class

#-------------------------------------------------------------------------------------------------#

@njit()
def Seeds_3D(Field: np.ndarray, dx: float, www: float, aaa: float):
    """
    Generates random 3D spherical seeds.

    Args:
        Field (np.ndarray): Array of phase fields.
        dx (float): Grid size.
        www (float): _description_
        aaa (float): _description_
    """
    
    NGrains, Nx, Ny, Nz = Field.shape
    
    r_nuclei = 3*dx
    
    for i in range(1, NGrains):   # loop through grains
        
        x_nuclei = int(np.random.rand()*Nx)
        y_nuclei = int(np.random.rand()*Ny)
        z_nuclei = int(np.random.rand()*Nz)
            
        for k in range(Nz):   # loop through z
            for m in range(Ny):   # loop through y
                for l in range(Nx):   # loop through x

                    # PBCs
                    if l > Nx-1: 
                        l = l - Nx
                    if l < 0:
                        l = l + Nx
                    if m > Ny-1:
                        m = m - Ny
                    if m < 0:
                        m = m + Ny
                    if k > Nz-1:
                        k = k - Nz
                    if k < 0:
                        k = k + Nz

                    r = np.sqrt((l*dx-x_nuclei*dx)**2 + \
                                (m*dx-y_nuclei*dx)**2 + \
                                (k*dx-z_nuclei*dx)**2) - r_nuclei

                    tmp = np.sqrt(2.*www)/aaa*r
                    phi_tmp = 0.5*(1.-np.sin(tmp))

                    if tmp >= np.pi/2.:
                        phi_tmp = 0.

                    if phi_tmp > 0: 
                        Field[i,l,m,k] = phi_tmp            
                        Field[0,l,m,k] = Field[0,l,m,k]-Field[i,l,m,k]  

                    # end l
                # end m
            # end k

    # end i

#-------------------------------------------------------------------------------------------------#
    
@njit()
def Seeds_2D(Field: np.ndarray, dx: float, www: float, aaa: float, nSeeds: int):
    """
    Generated random 2D circular seeds.

    Args:
        Field (np.ndarray): Array of phase fields.
        dx (float): Grid size.
        www (float): _description_
        aaa (float): _description_
    """
    
    NGrains, Nx, Ny = Field.shape
        
    r_nuclei = 3*dx
    
    for i in range(1, nSeeds):  # loop through grains
            
                x_nuclei = int(np.random.rand()*Nx)
                y_nuclei = int(np.random.rand()*Ny)
                    
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
                            Field[i,l,m] = phi_tmp            
                            Field[0,l,m] = Field[0,l,m] - Field[i,l,m]  

                    # end l
                        
                # end m

            # end i
            
#-------------------------------------------------------------------------------------------------#