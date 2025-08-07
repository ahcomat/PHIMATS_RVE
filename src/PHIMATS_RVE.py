#-----------------------------------------------------------------------------#
#																			  #
# Created    : 2022															  #
# Author(s)  : Abdelrahman Hussein a.h.a.hussein@outlook.com				  #
#																			  #
# This file is part of the software `PHIMATS RVE`							  #
#																			  #
#-----------------------------------------------------------------------------#

import sys
import os
import numpy as np
import h5py
from datetime import datetime

from numba_kernels import *

class PhaseField:
    """
    A class for managing phase-field evolution.
    """

    def __init__(self, InputData):
        
        self.Simul    = InputData["Simul"] 
        self.Nx       = InputData["Nx"]
        self.Ny       = InputData["Ny"]
        self.Nz       = InputData["Nz"]
        self.dx       = InputData["dx"]
        self.iWidth   = InputData["iWidth"]
        self.dt       = InputData["dt"]
        self.sigma    = InputData["sigma"]
        self.eee      = InputData["eee"]
        self.NumberOfGrains = InputData["NumberOfGrains"]

        self.eta   = self.iWidth*self.dx
        self.aaa   = 2.0 / np.pi * np.sqrt(2.0*self.eta*self.sigma)
        self.www   = 4.0 * self.sigma/self.eta
        self.pmobi = np.pi**2/(8.*self.eta)*3.0e-13
        
        # Check if 3D
        if self.Nz == 0:
            self.is3D = False
        else:
            self.is3D = True

        # Arrays for storing phase field data
        if not self.is3D:
            
            self.PhaseFields = np.zeros((self.NumberOfGrains, self.Nx, self.Ny))
            self.PhaseFieldsDot = np.zeros((self.NumberOfGrains, self.Nx, self.Ny))
            self.PhaseFieldID = np.zeros((self.Nx, self.Ny))
            self.gPhi = np.zeros((self.Nx, self.Ny))
            # Arrays for updating grain ID and number of grains 
            self.GrainID = np.zeros((self.NumberOfGrains, self.Nx, self.Ny), dtype=int)
            self.GrainNum = np.zeros((self.Nx, self.Ny), dtype=int)
            
        else :
            
            self.PhaseFields = np.zeros((self.NumberOfGrains, self.Nx, self.Ny, self.Nz))
            self.PhaseFieldsDot = np.zeros((self.NumberOfGrains, self.Nx, self.Ny, self.Nz))
            self.PhaseFieldID = np.zeros((self.Nx, self.Ny, self.Nz))
            self.gPhi = np.zeros((self.Nx, self.Ny, self.Nz))
            # Arrays for updating grain ID and number of grains 
            self.GrainID = np.zeros((self.NumberOfGrains, self.Nx, self.Ny, self.Nz), dtype=int)
            self.GrainNum = np.zeros((self.Nx, self.Ny, self.Nz), dtype=int)

        # Initialzie mobility matrix
        self.InitializeMatrices()
        
        # Write summary of simulation data ------
        print("\nSimulation: ", self.Simul)
        print("---------------------------------------")
        now = datetime.now()
        dt_string = now.strftime("%d-%m-%Y %H:%M:%S")
        print("Created on:", dt_string)
        print("Simulation size: ",  self.Nx, " x ", self.Ny, " x ", self.Nz, "grid points")
        print("Spacing dx: %.4e m" % self.dx)
        print("Numerical time increment dt: %.4e \n" % self.dt)
        
        # self.WriteSummary()

#-------------------------------------------------------------------------------------------------#
    
    def WriteSummary(self):
        """
        Writes simulation data summary to a txt file. 
        """
        
        original_stdout = sys.stdout
        
        with open(self.Simul+'_Data.txt', 'w') as f:
            sys.stdout = f
            print("Simulation: ", self.Simul)
            print()
            print("---------------------------------------")
            print()
            now = datetime.now()
            dt_string = now.strftime("%d-%m-%Y %H:%M:%S")
            print("Created on:", dt_string)
            print("Simulation size: ",  self.Nx, " x ", self.Ny, " x ", self.Nz, "grid points")
            print("Spacing dx: %.4e m" % self.dx)
            print("Numerical time increment dt: %.4e" % self.dt)
            sys.stdout = original_stdout
            
#-------------------------------------------------------------------------------------------------#      

    def InitializeMatrices(self):
    
        self.wij = np.zeros((self.NumberOfGrains, self.NumberOfGrains))
        self.aij = np.zeros((self.NumberOfGrains, self.NumberOfGrains))
        self.mij = np.zeros((self.NumberOfGrains, self.NumberOfGrains))
        self.eij = np.zeros((self.NumberOfGrains, self.NumberOfGrains))
        
        for i in range(self.NumberOfGrains):
            for j in range(self.NumberOfGrains):

                self.wij[i, j] = self.www
                self.aij[i, j] = self.aaa
                self.mij[i, j] = self.pmobi
                self.eij[i,j] = 0.0

                if i == j:
                    self.wij[i, j] = 0
                    self.aij[i, j] = 0
                    self.mij[i, j] = 0

                if i == 0 or j == 0:
                    self.eij[i,j] = self.eee
                if i < j:
                    self.eij[i,j] = -self.eij[i,j]

#-----------------------------------------------------------------------------#

    def OpenFileHDF5(self, overwrite=False):
        """
        Creates hdf5 file for storing simulation results. If the file exists, behavior depends on `overwrite`.

        Args:
            overwrite (bool): Overwrite the file if it exists. Defaults to True.

        Raises:
            OSError: If the file exists and overwrite is set to False.
        """
        file_path = self.Simul + ".hdf5"

        try:
            
            if overwrite:
                self.fh5 = h5py.File(file_path, "w")  # overwrite
            else:
                self.fh5 = h5py.File(file_path, "a")  # append if exists, or create new
                
        except OSError as e:
            raise RuntimeError(f"Failed to open HDF5 file '{file_path}': {e}")
        
#-----------------------------------------------------------------------------#

    def CloseFileHDF5(self):
        """
        Closes the currently open hdf5 file.

        Raises:
            RuntimeError: If no file is open or an error occurs during closing.
        """
        try:
            if hasattr(self, "fh5") and self.fh5:
                self.fh5.close()
                # print("HDF5 file closed successfully.")
            else:
                print("No HDF5 file is currently open.")
        except Exception as e:
            raise RuntimeError(f"Error occurred while closing the HDF5 file: {e}")
    
#-----------------------------------------------------------------------------#

    def WriteInputData(self):
        
        # Open hdf5 file 
        self.OpenFileHDF5()

        # SimulationParameters -------

        self.fh5.attrs["Simulation"] =  self.Simul
        
        self.grp_Sim_Params = self.fh5.create_group('SimulationParameters')
        
        self.grp_Sim_Params.create_dataset("Nx", data=self.Nx, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("Ny", data=self.Ny, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("Nz", data=self.Nz, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("NumberOfGrains", data=self.NumberOfGrains, dtype = np.int64)
        
        if self.is3D:
            self.grp_Sim_Params.create_dataset("is3D", data=1, dtype = np.int64)
        else:
            self.grp_Sim_Params.create_dataset("is3D", data=0, dtype = np.int64)

        self.grp_Sim_Params.create_dataset("dx", data=self.dx, dtype = np.float64)
        self.grp_Sim_Params.create_dataset("dt", data=self.dt, dtype = np.float64)

        self.grp_Sim_Params.create_dataset("wij", data=self.wij, dtype = np.float64, compression="gzip", compression_opts=5) 
        self.grp_Sim_Params.create_dataset("aij", data=self.aij, dtype = np.float64, compression="gzip", compression_opts=5) 
        self.grp_Sim_Params.create_dataset("mij", data=self.mij, dtype = np.float64, compression="gzip", compression_opts=5) 
        self.grp_Sim_Params.create_dataset("eij", data=self.eij, dtype = np.float64, compression="gzip", compression_opts=5) 
        
        self.grp_Arrays = self.fh5.create_group('Arrays')
        self.grp_Arrays.attrs["Contains"] = "Full-field arrays"
        
        self.grp_PhaseFields = self.grp_Arrays.create_group('PhaseFields')
        self.grp_PhaseFields.attrs["Contains"] = "Phase-field IDs time series"
        
        self.grp_gPhi = self.grp_Arrays.create_group('gPhi')
        self.grp_gPhi.attrs["Contains"] = "Interfaces gPhi time series"
        
        self.grp_phiID = self.grp_Arrays.create_group('phiID')
        self.grp_phiID.attrs["Contains"] = "Phase-field IDs time series"
        
        now = datetime.now()
        dt_string = now.strftime("%d-%m-%Y %H:%M:%S")
        
        if self.is3D:
            
            set_phiID = self.grp_phiID.create_dataset("phiID_"+str(0), data=self.PhaseFieldID, 
            dtype=np.float64, compression="gzip", compression_opts=5)
            set_phiID.attrs["Written on"] = dt_string
            set_gPhi= self.grp_gPhi.create_dataset("gPhi_"+str(0), data=self.gPhi, dtype=np.float64,
            compression="gzip", compression_opts=5)
            set_gPhi.attrs["Written on"] = dt_string
            
        else :
            
            set_phiID = self.grp_phiID.create_dataset("phiID_"+str(0), data=self.PhaseFieldID, 
            dtype=np.float64, compression="gzip", compression_opts=5)
            set_phiID.attrs["Written on"] = dt_string
            set_gPhi = self.grp_gPhi.create_dataset("gPhi_"+str(0), data=self.gPhi, dtype = np.float64,
            compression="gzip", compression_opts=5)
            set_gPhi.attrs["Written on"] = dt_string
            
            for i in range(self.NumberOfGrains):
                
                self.grp_PhaseFields.create_dataset("PhaseFields_"+str(i), data=self.PhaseFields[i,:,:], 
                dtype=np.float64, compression="gzip", compression_opts=5)
        
        self.CloseFileHDF5()
            
#-------------------------------------------------------------------------------------------------#
            
    def WriteResults(self, tStep):
        
        now = datetime.now()
        dt_string = now.strftime("%d-%m-%Y %H:%M:%S")
        
        self.OpenFileHDF5()
        
        if self.is3D:
            
            set_phiID = self.fh5.create_dataset("Arrays/phiID/phiID_"+str(tStep), data=self.PhaseFieldID, 
            dtype=np.float32, compression="gzip", 
            compression_opts=5)
            set_phiID.attrs["Written on"] = dt_string
            self.fh5.create_dataset("Arrays/gPhi/gPhi_"+str(tStep), data=self.gPhi, dtype=np.float32,
            compression="gzip", compression_opts=5)
            
        else :
            
            set_phiID = self.fh5.create_dataset("Arrays/phiID/phiID_"+str(tStep), data=self.PhaseFieldID, 
            dtype=np.float32, compression="gzip", compression_opts=5)
            set_phiID.attrs["Written on"] = dt_string
            set_gPhi = self.fh5.create_dataset("Arrays/gPhi/gPhi_"+str(tStep), data=self.gPhi, dtype = np.float32,
            compression="gzip", compression_opts=5)
            set_gPhi.attrs["Written on"] = dt_string
            
        self.CloseFileHDF5()
        
#-------------------------------------------------------------------------------------------------#

    def Calc_gPhi(self):
        """
        Calculates the diffuse interface function gPhi.
        """

        if self.is3D:
            self.gPhi = n_Calc_gPhi_3D(self.PhaseFields)
        else:
            self.gPhi = n_Calc_gPhi_2D(self.PhaseFields)
    
#-------------------------------------------------------------------------------------------------#

    def getPFID(self):
        """
        Evaluates the grain ID for the order parameters.
        """

        if self.is3D:
            
            n_getPFID_3D(self.PhaseFields, self.PhaseFieldID, self.GrainID, self.GrainNum)
        else:
            n_getPFID_2D(self.PhaseFields, self.PhaseFieldID, self.GrainID, self.GrainNum)

#-------------------------------------------------------------------------------------------------#

    def UpdateGrain(self):
        
        if self.is3D:
            n_UpdateGrain_3D(self.GrainID, self.GrainNum, self.PhaseFields)
        else:
            n_UpdateGrain_2D(self.GrainID, self.GrainNum, self.PhaseFields)
        
#-------------------------------------------------------------------------------------------------#

    def UpdateTimeStep(self):
        """
        Solve the phase-field evolution by explicit finite difference. 
        """

        if self.is3D:
            
            n_UpdateGrain_3D(self.GrainID, self.GrainNum, self.PhaseFields)
            n_CalcTimeStep_3D(self.PhaseFields, self.dx, self.dt, self.PhaseFieldsDot, self.GrainID, self.GrainNum, self.eij, self.wij, self.mij, self.aij)
        
        else :
            
            n_UpdateGrain_2D(self.GrainID, self.GrainNum, self.PhaseFields)
            n_CalcTimeStep_2D(self.PhaseFields, self.dx, self.dt, self.PhaseFieldsDot, self.GrainID, self.GrainNum, self.eij, self.wij, self.mij, self.aij)   

#-------------------------------------------------------------------------------------------------#