#!/bin/python3
# Import Libraries
import os
# import tempfile
from pylab import *

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *


# Set up the simulation
# Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WR90')
CURRENT_PATH = os.path.dirname(os.path.abspath(__file__))
APP_CSXCAD_PATH = '~/opt/openEMS/bin/'
Sim_Path = os.path.join(CURRENT_PATH, 'Rect_WR90')

post_proc_only = False
unit = 1e-3  # drawing unit in mm

# waveguide dimensions
# WR90
a = 23       # waveguide width
b = 10       # waveguide height
length = 50  # arms length

# frequency range of interest
f_start = 7e9
f_0     = 10e9
f_stop  = 12e9
lambda0 = C0/f_0/unit

# waveguide TE-mode definition
TE_mode = 'TE10'

# targeted mesh resolution
mesh_res = lambda0/30
print(f"{mesh_res=}")

# Setup FDTD parameter & excitation function
FDTD = openEMS(NrTS=1e5, EndCriteria=1e-4, TimeStepFactor=0.25)
FDTD.SetGaussExcite(0.5*(f_start+f_stop), 0.5*(f_stop-f_start))

# boundary conditions
FDTD.SetBoundaryCond([0, 0, 0, 0, 3, 3])

# Setup geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

x_start = 0
# x_end = a
# x_end = a + 1
x_end = a + length

y_start = 0
# y_end = b
# y_end = b + 1
y_end = b + length

z_start = 0 -length
z_end = b + length

mesh.AddLine('x', [x_start, x_end])
mesh.AddLine('y', [y_start, y_end])
mesh.AddLine('z', [z_start, z_end])

### Apply the waveguide port
### Port 1
ports = []
start = [0, 0,              z_start + 20]
stop  = [a, b, z_start + 20 + 5*mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(0, start, stop, 'z', a*unit, b*unit, "TE10", 1))

# Port 2
start = [0, 0,              z_end - 20]
stop  = [a, b, z_end - 20 - 5*mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(1, start, stop, 'z', a*unit, b*unit, "TE10"))

mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

### ADD PRIMITIVES ################################################################################
### Add waveguide walls
metal = CSX.AddMetal('metal') # create a metal property with name "metal"
start = [x_start, y_start, z_start]
stop  = [x_end,   y_end,     z_end]
box = metal.AddBox(start, stop)
box.SetPriority(1)

### Add waveguide space
air = CSX.AddMaterial('air', epsilon=1)
start = [x_start, y_start, z_start]
stop  = [a,       b,         z_end]
box = air.AddBox(start, stop)
box.SetPriority(2)

# Smooth mesh
FDTD.AddEdges2Grid(dirs='all', properties=air, primitives=box)
mesh.AddLine('x', [23 - 0.1, 23 + 0.1])
mesh.AddLine('y', [10 - 0.1, 10 + 0.1])
mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

# Run the simulation ##############################################################################
if True:  # debugging only
    CSX_file = os.path.join(Sim_Path, 'rect_wg.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(APP_CSXCAD_PATH + AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

# Postprocessing & plotting
freq = linspace(f_start, f_stop, 201)
for port in ports:
    port.CalcPort(Sim_Path, freq)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc

print(f"{abs(s11)=}")
print(f"{abs(s21)=}")

# Plot s-parameter
plt.figure()
plt.subplot(1, 2, 1)
plt.plot(freq*1e-6, 20*log10(abs(s11)), 'r--', marker="1", markevery=50, linewidth=2, label='$S_{11}$')
plt.grid()
plt.plot(freq*1e-6, 20*log10(abs(s21)), 'g--', marker="2", linewidth=2, markevery=55, label='$S_{21}$')
plt.legend()
plt.ylabel('S-Parameter (dB)')
plt.xlabel(r'frequency (MHz) $\rightarrow$')

plt.subplot(1, 2, 2)
plt.plot(freq*1e-6, abs(s11**2), 'r--', linewidth=2, label='$S_{11}^{2}$')
plt.grid()
plt.plot(freq*1e-6, abs(s21**2), 'g--', linewidth=2, label='$S_{21}^{2}$')
plt.plot(freq*1e-6, abs(s11)**2 + abs(s21)**2, 'k-', linewidth=2, label='$Summ$')
plt.legend()
plt.ylabel('S-Parameter')
plt.xlabel(r'frequency (MHz) $\rightarrow$')
plt.show()
