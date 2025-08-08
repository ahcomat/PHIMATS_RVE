import numpy as np
import h5py
import gmsh

from scipy.interpolate import RectBivariateSpline

#-----------------------------------------------------------------------------#

def Resample2D(Cx: np.ndarray, Mx: int, dx: float = 0.25e-6) -> np.ndarray :
    """
    Resamples a 2D array

    Args:
        Cx (np.ndarray): Array to be interpolated (Should be square)
        Mx (int): New array dimension
        dx (float, optional): Grid size. Defaults to 0.25e-6.

    Raises:
        ValueError: If Cx is not a square array. 

    Returns:
        np.ndarray: Interpolated Cx to new size [Mx, Mx]
    """
    
    if Cx.shape[0] != Cx.shape[1]:
        raise ValueError(f"Expected a square matrix, but got shape {Cx.shape}")
    
    x = np.linspace(0, Cx.shape[0]*dx, Cx.shape[0])

    interp_spline2 = RectBivariateSpline(x, x, Cx)

    xx = np.linspace(0, x.max(), Mx)
    Cx_interp = interp_spline2(xx,xx)
    Cx_interp[Cx_interp<=0] = 0 # Make sure no negative values

    return Cx_interp

# -------------------------------------------------------------------------------------------------#

def CreateStructuredMesh(Lx: float, Mx: int, filename: str = "RVE_Structured", dim: int = 2, export_format: str = "inp"):
	"""
	Creates a structured mesh for the RVE using Gmsh API.
 
	Args:
		Lx (float): Edge length of the RVE [m]
		Mx (int): Number of grid pints
		filename (str, optional): Resulting file name. Defaults to "RVE".
		dim (int, optional): Dimension. Defaults to 2.
		export_format (str, optional): File format for the exported mesh. **Note:** PHIMATS currently supports only the .inp extension. Defaults to"inp".
	"""
	
	gmsh.initialize()
	gmsh.model.add("RVE_Model")

	# Define points
	p1 = gmsh.model.geo.addPoint(Lx,   0, 0, 1.0)
	p2 = gmsh.model.geo.addPoint(Lx, Lx, 0, 1.0)
	p3 = gmsh.model.geo.addPoint( 0, Lx, 0, 1.0)
	p4 = gmsh.model.geo.addPoint( 0,   0, 0, 1.0)

	# Define lines
	l1 = gmsh.model.geo.addLine(p4, p1)
	l2 = gmsh.model.geo.addLine(p1, p2)
	l3 = gmsh.model.geo.addLine(p2, p3)
	l4 = gmsh.model.geo.addLine(p3, p4)

	# Curve loop and surface
	loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
	surface = gmsh.model.geo.addPlaneSurface([loop])

	# Transfinite definition (structured grid)
	gmsh.model.geo.mesh.setTransfiniteCurve(l1, Mx)
	gmsh.model.geo.mesh.setTransfiniteCurve(l2, Mx)
	gmsh.model.geo.mesh.setTransfiniteCurve(l3, Mx)
	gmsh.model.geo.mesh.setTransfiniteCurve(l4, Mx)
	gmsh.model.geo.mesh.setTransfiniteSurface(surface)

	gmsh.model.geo.synchronize()

	# Physical group
	gmsh.model.addPhysicalGroup(2, [surface], tag=1)
	gmsh.model.setPhysicalName(2, 1, "RVE")

	gmsh.model.mesh.generate(dim)

	gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)

	# Export
	gmsh.write(f"{filename}.{export_format}")

	gmsh.finalize()

	print(f"Mesh saved to '{filename}.{export_format}'")

# -------------------------------------------------------------------------------------------------#

def BackgroundMesh(Lx, pos_file="nodal_values.pos", msh_file="RVE_Adaptive.inp"):
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("RVE_geo")

    # Geometry 
    p1 = gmsh.model.geo.addPoint(Lx,   0, 0, 1.0)
    p2 = gmsh.model.geo.addPoint(Lx, Lx, 0, 1.0)
    p3 = gmsh.model.geo.addPoint( 0, Lx, 0, 1.0)
    p4 = gmsh.model.geo.addPoint( 0,   0, 0, 1.0)

    l1 = gmsh.model.geo.addLine(p4, p1)
    l2 = gmsh.model.geo.addLine(p1, p2)
    l3 = gmsh.model.geo.addLine(p2, p3)
    l4 = gmsh.model.geo.addLine(p3, p4)

    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surface = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.synchronize()
    
	# Physical group
    gmsh.model.addPhysicalGroup(2, [surface], tag=1)
    gmsh.model.setPhysicalName(2, 1, "RVE")

	# Merge background mesh from .pos file 
    gmsh.merge(pos_file)
    bg_field = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setNumber(bg_field, "ViewIndex", 0)
    gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

	# Mesh options
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 1)  # MeshAdapt
    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1) # Export groups of nodes

    # Generate and export 
    gmsh.model.mesh.generate(2)
    gmsh.write(msh_file)
    gmsh.finalize()

    print(f"Mesh generated and saved to: {msh_file}")

# -------------------------------------------------------------------------------------------------#
