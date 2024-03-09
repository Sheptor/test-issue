# Import Libraries
import os
# import tempfile
from pylab import *

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
# import CSXCAD

# Set up the simulation
# Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WR90')
CURRENT_PATH = os.path.dirname(os.path.abspath(__file__))
APP_CSXCAD_PATH = '/home/semicond/opt/openEMS/bin/'
Sim_Path = os.path.join(CURRENT_PATH, 'Rect_WR90')

post_proc_only = False
unit = 1e-3  # drawing unit in um

# waveguide dimensions
# WR90
a = 23        # waveguide width
b = 10        # waveguide height
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
# FDTD.SetBoundaryCond([0, 0, 0, 0, 3, 3])
FDTD.SetBoundaryCond(["PEC", "PML_8", "PEC", "PML_8", "PML_8", "PML_8"])

# Setup geometry & mesh
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

x_start = 0
x_end = a + length

y_start = 0
y_end = b + length

z_start = 0 - length
z_end = b + length

mesh.AddLine('x', [x_start, x_end])
mesh.AddLine('y', [y_start, y_end])
mesh.AddLine('z', [z_start, z_end])

#  Apply the waveguide port

#  Port 1
ports = []
start = [0, 0,              z_start + 20]
stop  = [a, b, z_start + 20 + 5*mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(0, start, stop, 'z', a*unit, b*unit, "TE10", 1))

# Port 2
start = [0, 0,                z_end - 20]
stop  = [a, b, z_end - 20 - 5 * mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(1, start, stop, 'z', a*unit, b*unit, "TE10"))

# Port 3 (E-Arm)
start = [0, y_end - 20 - 5 * mesh_res, -b/2]
stop  = [a, y_end - 20,                 b/2]
mesh.AddLine('y', [start[1], stop[1]])
ports.append(FDTD.AddRectWaveGuidePort(2, start, stop, 'y', a*unit, b*unit, "TE10"))

# Port 2 (H-Arm)
start = [x_end - 20,                0, -a/2]
stop  = [x_end - 20 - 5 * mesh_res, b,  a/2]
mesh.AddLine('x', [start[0], stop[0]])
ports.append(FDTD.AddRectWaveGuidePort(3, start, stop, 'x', a*unit, b*unit, "TE10"))

mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

# ADD PRIMITIVES #######################################################################################################

# Add waveguide walls
metal = CSX.AddMetal('metal')  # create a metal property with name "metal"
start = [x_start, y_start, z_start]
stop  = [x_end,   y_end,     z_end]
box = metal.AddBox(start, stop)
box.SetPriority(1)

# Add waveguide space
co_linear_air = CSX.AddMaterial('co_linear_air', epsilon=1)
start = [x_start, y_start, z_start]
stop  = [a,       b,         z_end]
co_linear = co_linear_air.AddBox(start, stop)
co_linear.SetPriority(2)


e_air = CSX.AddMaterial('e_air', epsilon=1)
start = [0, 0,     -b/2]
stop  = [a, y_end,  b/2]
e_arm = e_air.AddBox(start, stop)
e_arm.SetPriority(3)

h_air = CSX.AddMaterial('h_air', epsilon=1)
start = [0,     0, -a/2]
stop  = [x_end, b,  a/2]
h_arm = h_air.AddBox(start, stop)
h_arm.SetPriority(4)

# Run the simulation ###################################################################################################
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
freq = linspace(f_start, f_stop, 2001)
for port in ports:
    port.CalcPort(Sim_Path, freq)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc
s31 = ports[2].uf_ref / ports[0].uf_inc
s41 = ports[3].uf_ref / ports[0].uf_inc
summ = abs(s11)**2 + abs(s21)**2 + abs(s31)**2 + abs(s41)**2

print(f"{abs(s11)=}")
print(f"{abs(s21)=}")
print(f"{abs(s31)=}")
print(f"{abs(s41)=}")
print(f"{summ=}")

# Plot s-parameter

plt.figure()
plt.subplot(1, 2, 1)
plt.plot(freq*1e-6, 20*log10(abs(s11)), 'r--', marker="1", markevery=50, linewidth=2, label='$S_{11}$')
plt.grid()
plt.plot(freq*1e-6, 20*log10(abs(s21)), 'g--', marker="2", linewidth=2, markevery=55, label='$S_{21}$')
plt.plot(freq*1e-6, 20*log10(abs(s31)), 'b--', marker="2", linewidth=2, markevery=55, label='$S_{31}$')
plt.plot(freq*1e-6, 20*log10(abs(s41)), 'y--', marker="2", linewidth=2, markevery=55, label='$S_{41}$')
plt.legend()
plt.ylabel('S-Parameter (dB)')
plt.xlabel(r'frequency (MHz) $\rightarrow$')

plt.subplot(1, 2, 2)
plt.plot(freq*1e-6, abs(s11**2), 'r--', linewidth=2, label='$S_{11}^{2}$')
plt.grid()
plt.plot(freq*1e-6, abs(s21**2), 'g--', linewidth=2, label='$S_{21}^{2}$')
plt.plot(freq*1e-6, abs(s31**2), 'b--', linewidth=2, label='$S_{31}^{2}$')
plt.plot(freq*1e-6, abs(s41**2), 'y--', linewidth=2, label='$S_{41}^{2}$')
plt.plot(freq*1e-6, summ, 'k-', linewidth=1, label='$Summ$')
plt.legend()
plt.ylabel('S-Parameter')
plt.xlabel(r'frequency (MHz) $\rightarrow$')
plt.show()
