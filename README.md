# CALiPPSO.py - Jamming Hard Spheres in Python
This is a draft of the python version of "CALiPPSO.jl: A Linear Programming Algorithm for Jamming Hard Spheres".
</br>
At the moment, the 2d version has been tested(while still not optimally working), while the 3d version is usable only up to the creation of the pre-jammed enviroment.

## Installation
The code is written in Python 3.11.4 and requires the following packages:
- numpy
- scipy
- matplotlib
- gurobipy

To install the packages, run the following command in the terminal:
You can install these packages using `pip` with the following command:

```
pip install numpy scipy matplotlib gurobipy
```

## Usage
The code is run by executing the file `main.py` in the terminal.

Run the code with:

```
python3 src/main.py
```

followed by these possible options:

* `--num_spheres`: Number of spheres to generate and jam (default=5)
* `--plot`: Whether to plot the initial and final configurations (default=1, 0 to disable)

