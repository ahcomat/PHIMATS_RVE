<img src="PHIMATS_RVE.png" width="200">

PHIMATS *RVE* is a package for generating representative volume elements (RVEs) based on the multiphase-field method. While it can be used for general phase transformation based on grain growth, it is primarily developed to generate RVEs for the FEM code **[PHIMATS](https://github.com/ahcomat/PHIMATS.git)**.

## Installation Guide

### **1️⃣ System Requirements**
✅ **Linux OS required** (Recommended: Ubuntu, Debian, Fedora)
✅ **Windows:** Use **Windows Subsystem for Linux (WSL)**
⚠️ **macOS:** Potentially supported but untested 

---

### **2️⃣ Python Dependencies (via Conda)**

All required Python libraries are included in the **`phimats_rve`** Conda environment.

To create the environment, navigate to `PHIMATS_RVE` directory and run:
```sh
conda env create -f environment.yml
```

---

### **3️⃣ Setting Up Environment Variables**
Add the following lines **(after modifying the paths)** to your `~/.bashrc` file:
```sh

export PYTHONPATH=$PYTHONPATH:/path_to_PHIMATS_RVE/src/

export HDF5_USE_FILE_LOCKING=FALSE
```
After adding these, apply the changes by running:
```sh
source ~/.bashrc
```

--- 

### Acknowledgment

This implementation is partially based on the original multiphase-field code published by the [Yamanaka Research Group, Tokyo University of Agriculture and Technology](https://web.tuat.ac.jp/~yamanaka/opensource.html). The original code served as a starting point for developing the MPF functionality in PHIMATS *RVE*. All rights to the original code remain with its respective developers, and use is subject to their [terms of use](https://web.tuat.ac.jp/~yamanaka/opensource.html).

---

### License
PHIMATS *RVE* is released under the **[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html)**, or later.  

