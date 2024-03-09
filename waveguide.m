#!/bin/octave


close all
clear
clc

%% Setup the Simulation
physical_constants;
unit = 1e-3; % all length in mm
openEMS_opts = '';
openEMS_opts = [openEMS_opts ' --disable-dumps'];
openEMS_opts = [openEMS_opts ' --dump-statistics'];

% openEMS_opts = [openEMS_opts ' --debug-material'];
% openEMS_opts = [openEMS_opts ' --engine=basic'];

% openEMS_opts = [openEMS_opts ' --debug-PEC'];
% openEMS_opts = [openEMS_opts ' --debug-operator'];
% openEMS_opts = [openEMS_opts ' --debug-boxes'];
% openEMS_opts = [openEMS_opts ' --debug-CSX'];
Settings = [];
Settings.LogFile = 'openEMS.log';


a = 23.0;
b = 10.0;
c = 60.0;

f_start = 7e9;
f_stop  = 13e9;

FDTD = InitFDTD('NrTS',3e5, 'EndCriteria', 1e-5, 'OverSampling', 5,  'TimeStepFactor', 0.4);
FDTD = SetGaussExcite(FDTD,0.5*(f_start+f_stop),0.5*(f_stop-f_start));
BC = {'PEC' 'MUR' 'PEC' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond(FDTD, BC);


CSX = InitCSX();
mesh_res=0.5
mesh.x = [0:mesh_res:c];
mesh.y = [0:mesh_res:c];
mesh.z = [-c:mesh_res:c];

CSX = DefineRectGrid(CSX, unit, mesh);

CSX = AddMetal(CSX, 'space');
start = [mesh.x(1) mesh.y(1) mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddBox(CSX,'space',0,start,stop);

% Create Substrate
CSX = AddMaterial(CSX, 'vacuum');
CSX = SetMaterialProperty( CSX, 'vacuum', 'Epsilon', 1);
start = [0 0 mesh.z(1)];
stop  = [a b mesh.z(end)];
CSX = AddBox(CSX, 'vacuum', 1, start, stop);


start = [mesh.x(1) 0 -a/2];
stop  = [mesh.x(end) b a/2];
CSX = AddBox(CSX, 'vacuum', 1, start, stop);


start = [0 mesh.y(1) -b/2];
stop  = [a mesh.y(end) b/2];
CSX = AddBox(CSX, 'vacuum', 1, start, stop);


%% apply the waveguide port %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
port_z = 15;
length_port = int64(port_z/(mesh_res));

index_a = int64(a/(mesh_res));
index_b = int64(b/(mesh_res));
index_c = int64(c/(mesh_res));
%waveguide TE-mode definition
TE_mode = 'TE10';
ext_port = [1 0 0 0];

start=[mesh.x(1)   mesh.y(1)   mesh.z(length_port)];
stop =[mesh.x(index_a+1) mesh.y(index_b+1) mesh.z(length_port+10)];

start=[0  0  mesh.z(length_port)];
stop =[23 10 mesh.z(length_port+10)];

[CSX, port{1}] = AddRectWaveGuidePort( CSX, 0, 1, start, stop, 'z', a*unit, b*unit, TE_mode, ext_port(1));

start=[mesh.x(1)   mesh.y(1)   mesh.z(end-length_port)];
stop =[mesh.x(index_a+1) mesh.y(index_b+1) mesh.z(end-length_port-10)];

start=[0   0   mesh.z(end-length_port)];
stop =[23 10 mesh.z(end-length_port-10)];
[CSX, port{2}] = AddRectWaveGuidePort( CSX, 0, 2, start, stop, 'z', a*unit, b*unit, TE_mode, ext_port(2));

start=[mesh.x(end-length_port)   mesh.y(1)   mesh.z(index_c-index_a/2+1)];
stop =[mesh.x(end-length_port-10) mesh.y(index_b+1) mesh.z(index_c+index_a/2+1)];
[CSX, port{3}] = AddRectWaveGuidePort( CSX, 0, 3, start, stop, 'x', a*unit, b*unit, TE_mode, ext_port(3));

start=[mesh.x(1)    mesh.y(end-length_port)    mesh.z(index_c-index_b/2+1)];
stop =[mesh.x(index_a+1) mesh.y(end-length_port-10)   mesh.z(index_c+index_b/2+1)];
[CSX, port{4}] = AddRectWaveGuidePort( CSX, 0, 4, start, stop, 'y', a*unit, b*unit, TE_mode, ext_port(4));



%% Prepare and Run Simulation
Sim_Path = 'tmp_mod';
Sim_CSX = 'patch_ant.xml';

% create an empty working directory
[status, message, messageid] = rmdir(Sim_Path, 's'); % clear previous directory
[status, message, messageid] = mkdir(Sim_Path); % create empty simulation folder

% write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

% show the structure
CSXGeomPlot([Sim_Path '/' Sim_CSX]);
% exit()
RunOpenEMS(Sim_Path, Sim_CSX, [openEMS_opts ' '], Settings)



%% postproc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = linspace(f_start,f_stop,2000);
port = calcPort(port, Sim_Path, freq);

% s = [];
% for i = 1:1:size(ext_port)(0)
%
s11 = port{1}.uf.ref./ port{1}.uf.inc;
s21 = port{2}.uf.ref./ port{1}.uf.inc;
s31 = port{3}.uf.ref./ port{1}.uf.inc;
s41 = port{4}.uf.ref./ port{1}.uf.inc;
ZL = port{1}.uf.tot./port{1}.if.tot;
ZL_a = port{1}.ZL; % analytic waveguide impedance

name_obj = ['123' '123']
images = 'images';
[status, message, messageid] = mkdir(images);
%% plot s-parameter dB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hf = figure('visible','off');
plot(freq*1e-6,20*log10(abs(s11)),'k-','Linewidth',4);
xlim([freq(1) freq(end)]*1e-6);
grid on;
hold on;
plot(freq*1e-6,20*log10(abs(s21)),'r--','Linewidth',4);
plot(freq*1e-6,20*log10(abs(s31)),'g--','Linewidth',4);
plot(freq*1e-6,20*log10(abs(s41)),'b--','Linewidth',4);
l = legend('S_{11}','S_{21}','S_{31}','S_{41}','Location','northeast');
set(l,'FontSize',8);
ylabel('S-Parameter (dB)','FontSize',8);
xlabel('frequency (MHz) \rightarrow','FontSize',8);
saveas(hf, [images '/' 'S-param_dB_' name_obj '.png'], 'png');
delete(hf);

%% plot s-parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hf = figure('visible','off');
plot(freq*1e-6,abs(s11),'k-','Linewidth',4);
xlim([freq(1) freq(end)]*1e-6);
grid on;
hold on;
plot(freq*1e-6,abs(s21),'r--','Linewidth',4);
plot(freq*1e-6,abs(s31),'g--','Linewidth',4);
plot(freq*1e-6,abs(s41),'b--','Linewidth',4);
l = legend('S_{11}','S_{21}','S_{31}','S_{41}','Location','northeast');
set(l,'FontSize',8);
ylabel('S-Parameter','FontSize',8);
xlabel('frequency (MHz) \rightarrow','FontSize',8);
saveas(hf, [images '/' 'S-param_' name_obj  '.png'], 'png');
delete(hf);
