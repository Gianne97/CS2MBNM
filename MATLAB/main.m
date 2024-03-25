% Giacomo Giannetti, University of Florence
% 2024, January, 3
% Comparative study of two multipole-based numerical methods for 3D field-translation schemes

% Global system: the one used in CST
% Antenna system: the one for the sphere enclosing the antenna
% Scatterer system: the one of the scatterer

clear
close all
clc

PathKiel = "C:\Users\gigi\";
PathHome = "D:\";
PathUnifi= "C:\Users\giannetti\";
RootPath = PathUnifi;
% addpath(strcat(RootPath, "OneDrive - unifi.it\Università\Software\MATLAB\UsefulFunctions"))
addpath(strcat(RootPath, "OneDrive - unifi.it\Università\Dottorato\Lavori\NACS\MATLAB\3D\github_repo"))

ScatterSelection = 1; % 1 simple sphere, 2 multilayered sphere
SimSelector = 2; % 1, only antenna; 2, no scatterer; 3, with scatterer
Software = "CST"; % CST or HFSS

c = physconst('LightSpeed'); SaC.c = c; % Speed of light [m/s]
mu0 = 4*pi*1e-7; SaC.mu0 = mu0; % vacuum permeability [H/m]
eps0 = 1/(mu0*c^2); SaC.eps0 = eps0; % vacuum permittivity [F/m]
Z0 = sqrt(mu0/eps0); SaC.Z0 = Z0; % vacuum impedance [ohm]
Y0 = 1/Z0; SaC.Y0 = Y0; % vacuum admittance [ohm]

f0 = 2e9; SaC.f0 = f0;     % Working frequency [Hz]
w0 = 2*pi*f0; SaC.w0 = w0; % Working radiant frequency [rad/s]
k0 = w0/c; SaC.k0 = k0;     % Wavenumber in free space [1/m]
lambda = c/f0; % Wavelength [m]
typ = "h2"; SaC.typ = typ; % Type of spherical Hankel function

m2mm = 1e3;
mm2m = 1e-3;

% Scatterer center in the antenna reference system - case A
x0 = 15*mm2m; % [m]
y0 = 30*mm2m; % [m]
z0 =700*mm2m; % [m]
% z0 =900*mm2m; % [m]
C0 = [x0, y0, z0];
Rs = 30*mm2m; % sphere radius [m]
epsR = 2.2; % relative dielectric permittivity
[th0, ph0, r0] = CoorTran(C0(:));

% Antenna parameters
wg_width = 109.22*mm2m; % waveguide width, (m)
thickness = 2*mm2m; % wall thickness, (m)
wg_height = 54.61*mm2m; % waveguide height, (m)
horn_length = 575*mm2m; % horn length, (m)
wg_length = 75*mm2m; % waveguide length, (m)
taper_angle = 8.5; % (deg)

dR = 60*mm2m; % indicative extra diameter for the sphere enclosing the antenna (m)
drect = 20*mm2m; % distance from the horn for the rectangle enclosing the antenna (m)

a = 0.5*sqrt((x0 + (wg_width/2 + thickness))^2 + (y0 + (wg_height/2 + thickness))^2 + (z0 + (horn_length/2 + wg_length/2 + dR/2))^2); % radius of the minimum sphere centered in the scatterer reference frame and enclosing both objects (m)
Lmax1= ceil(k0*a + 7*(k0*a)^(1/3) + 3);
Lmax2= ceil(k0*a + 10*log(k0*a + pi));
ceil(k0*a);
ceil(k0*a + 10);
ceil(k0*a + 3*(k0*a)^(1/3));

xmin = -(0.5*wg_width +thickness+tan(deg2rad(taper_angle))*horn_length + drect);
xmax = -xmin;
ymin = -(0.5*wg_height+thickness+tan(deg2rad(taper_angle))*horn_length + drect);
ymax = -ymin;
zmin = -(thickness + drect);
zmax = wg_length+horn_length + drect;
SaC.xmin = xmin; SaC.xmax = xmax;
SaC.ymin = ymin; SaC.ymax = ymax;
SaC.zmin = zmin; SaC.zmax = zmax;

%% Pointlist for field extraction - Antenna Sphere with and without scatterer

N = 33; % "^m + 1
el = linspace(-pi/2, pi/2, N); % Elevation (rad)
az = linspace(-pi, pi, 2*N-1); % Azimuth (rad)
th = pi/2 - el; % Check MATLAB documentation
ph = az;
Ra = 385*mm2m; % Radius of the sphere enclosing the antenna (m)
zc = 355*mm2m; % Center of the sphere enclosing the antenna wrt global coordinate system (mm)

dd_el_az = diff(el(1:2));

% z0S= z0 - zc; % Center of the scatterer sphere wrt enclosing sphere
z0S=z0; % Center of the scatterer sphere wrt enclosing sphere

if Ra > norm([x0,y0,z0S]) - Rs
   warning("The scattering sphere intersects the sphere enclosing the antenna")
end

% [El, Az] = meshgrid(el, az);
% [Th, Ph] = meshgrid(th, ph);
% 
% Center = [0, 0, zc];
% PointList(Az, El, Ra, Center, m2mm, N, "PointlistAntennaSphere.txt")
% Center = [x0, y0, z0 + zc];
% PointList(Az, El, Rs, Center, m2mm, N, "PointlistScatterSphere.txt")

%% Pointlist antenna extraction GL

% I follow the indications I found in a paper (what paper?)
L = 60; % degree of the multipole expansion
QH = 'GL'; % Choose between different algorithms for the points of which the fields are extracted; 'GL' is for evenly-spaced points along the angular coordinates

Center = [0, 0, zc];
[N_ph, N_th] = PointList2(Ra, Center, m2mm, L, QH, "PointlistAntennaSphereGL.txt");

Center = [x0, y0, z0 + zc];
[~, ~] = PointList2(Rs, Center, m2mm, L, QH, "PointlistScatterSphereGL.txt");

%% Pointlist for field extraction - Antenna Rectangle with and without scatterer

WavelenRatio = 20;
dd = lambda/WavelenRatio;
Lx = xmax-xmin;
Ly = ymax-ymin;
Lz = zmax-zmin;

Nx = ceil(Lx/dd +1); Ny = ceil(Ly/dd +1); Nz = ceil(Lz/dd +1);
xR = linspace(xmin, xmax, Nx); dx = Lx/(Nx-1);
yR = linspace(ymin, ymax, Ny); dy = Ly/(Ny-1);
zR = linspace(zmin, zmax, Nz); dz = Lz/(Nz-1);

if not((dd > dx) && (dd > dy) && (dd > dz))
    error("There is an error")
end

[XR, YR, ZR] = meshgrid(xR, yR, zR);
XRG = reshape(XR, [],1);
YRG = reshape(YR, [],1);
ZRG = reshape(ZR, [],1);

% YZ planes
idxSides{1} = XRG == xmin;
idxSides{2} = XRG == xmax;
S{1} = [XRG(idxSides{1}), YRG(idxSides{1}), ZRG(idxSides{1})]; 
S{2} = [XRG(idxSides{2}), YRG(idxSides{2}), ZRG(idxSides{2})];
n{1} = "-x";
n{2} = "+x";

% XZ planes
idxSides{1} = YRG == ymin;
idxSides{2} = YRG == ymax;
S{3} = [XRG(idxSides{1}), YRG(idxSides{1}), ZRG(idxSides{1})];
S{4} = [XRG(idxSides{2}), YRG(idxSides{2}), ZRG(idxSides{2})];
n{3} = "-y";
n{4} = "+y";

% XY planes
idxSides{1} = ZRG == zmin;
idxSides{2} = ZRG == zmax;
S{5} = [XRG(idxSides{1}), YRG(idxSides{1}), ZRG(idxSides{1})];
S{6} = [XRG(idxSides{2}), YRG(idxSides{2}), ZRG(idxSides{2})];
n{5} = "-z";
n{6} = "+z";

Sides = [S{1}; S{2}; S{3}; S{4}; S{5}; S{6}]; % 2*(Nx*Ny + Ny*Nz + Ny*Nz)
SidesMin = [S{1}; S{3}; S{5}];
SidesMax = [S{2}; S{4}; S{6}];

SidesOrdered = m2mm*unique(Sides, 'rows', 'stable'); % 2*(Nx*Ny + Ny*Nz + Ny*Nz) - (2*2*(Nx + Ny) + 4*Nz) + 8

NpointsRect = size(SidesOrdered,1);

fid = fopen("FieldData\PointlistAntennaRectangle.txt", 'w+');
fprintf(fid, '%8.3f %8.3f %8.3f\n', transpose(SidesOrdered));
fclose(fid);

%% Field data

% Simulated fields on the enclosing box
ERectScatt0 = readtable("FieldData\Efield_2GHz_Rect_Scatterer0.txt");
HRectScatt0 = readtable("FieldData\Hfield_2GHz_Rect_Scatterer0.txt");
ERectScatt1 = readtable("FieldData\Efield_2GHz_Rect_Scatterer1.txt");
HRectScatt1 = readtable("FieldData\Hfield_2GHz_Rect_Scatterer1.txt");
ERS0 = table2array(ERectScatt0);
HRS0 = table2array(HRectScatt0);
ERS1 = table2array(ERectScatt1);
HRS1 = table2array(HRectScatt1);

% ESphereScatt0 = readtable("FieldData\Efield_2GHz_Sphere_NU_Scatterer0.txt");
% ESphereScatt0 = readtable("FieldData\Efield_2GHz_Sphere_Scatterer0.txt");
% HSphereScatt0 = readtable("FieldData\Hfield_2GHz_Sphere_Scatterer0.txt");
% ESphereScatt1 = readtable("FieldData\Efield_2GHz_Sphere_Scatterer1.txt");
% HSphereScatt1 = readtable("FieldData\Hfield_2GHz_Sphere_Scatterer1.txt");

% Simulated fields on the enclosing sphere
ESphereScatt0 = readtable("FieldData\Efield_2GHz_Sphere_GL_L60_Scatterer0.txt");
HSphereScatt0 = readtable("FieldData\Hfield_2GHz_Sphere_GL_L60_Scatterer0.txt");
ESphereScatt1 = readtable("FieldData\Efield_2GHz_Sphere_GL_L60_Scatterer1.txt");
HSphereScatt1 = readtable("FieldData\Hfield_2GHz_Sphere_GL_L60_Scatterer1.txt");
ESS0 = table2array(ESphereScatt0);
HSS0 = table2array(HSphereScatt0);
ESS1 = table2array(ESphereScatt1);
HSS1 = table2array(HSphereScatt1);

% Simulated fields on the scatterer surface
ESphereScattScatt0 = readtable("FieldData\Efield_2GHz_SphereScatt_GL_L60_Scatterer0.txt");
HSphereScattScatt0 = readtable("FieldData\Hfield_2GHz_SphereScatt_GL_L60_Scatterer0.txt");
ESSS0 = table2array(ESphereScattScatt0);
HSSS0 = table2array(HSphereScattScatt0);

% Field data extraction
time = 0.0e-9;
[Fields_RS0] = FieldsFromFWS(ERS0, HRS0, SaC);
Fields_RS0_t = TimeEvaluation(Fields_RS0, f0, time);
[Fields_RS1] = FieldsFromFWS(ERS1, HRS1, SaC);

Fields_RS1_t = TimeEvaluation(Fields_RS1, f0, time);
[Fields_SS0] = FieldsFromFWS(ESS0, HSS0, SaC);
Fields_SS0_t = TimeEvaluation(Fields_SS0, f0, time);

[Fields_SS1] = FieldsFromFWS(ESS1, HSS1, SaC);
Fields_SS1_t = TimeEvaluation(Fields_SS1, f0, time);

[Fields_SSS0] = FieldsFromFWS(ESSS0, HSSS0, SaC);
Fields_SSS0_t = TimeEvaluation(Fields_SSS0, f0, time);

xminSim = min(ERS0(:,1)); SaC.xminSim = xminSim;
xmaxSim = max(ERS0(:,1)); SaC.xmaxSim = xmaxSim;
yminSim = min(ERS0(:,2)); SaC.yminSim = yminSim;
ymaxSim = max(ERS0(:,2)); SaC.ymaxSim = ymaxSim;
zminSim = min(ERS0(:,3)); SaC.zminSim = zminSim;
zmaxSim = max(ERS0(:,3)); SaC.zmaxSim = zmaxSim;

% Equivalent currents
% SideStruct.
Sidx{1} = ERS0(:,1) == xminSim;
Slabel(1) = "xmin";
Sdd(1) = dy*dz;

Sidx{2} = ERS0(:,1) == xmaxSim;
Slabel(2) = "xmax";
Sdd(2) = dy*dz;

Sidx{3} = ERS0(:,2) == yminSim;
Slabel(3) = "ymin";
Sdd(3) = dx*dz;

Sidx{4} = ERS0(:,2) == ymaxSim;
Slabel(4) = "ymax";
Sdd(4) = dx*dz;

Sidx{5} = ERS0(:,3) == zminSim;
Slabel(5) = "zmin";
Sdd(5) = dx*dy;

Sidx{6} = ERS0(:,3) == zmaxSim;
Slabel(6) = "zmax";
Sdd(6) = dx*dy;

% Equivalent currents
[Je, Jm, JeR, JmR, xSides, ySides, zSidesGRF] = EquivalentCurrents(ERS0, Fields_RS0, Sidx);
zSides = zSidesGRF - zc*m2mm;

%% From Cartesian to spherical components - simulated electric field on the sphere

% output: reference frame of the sphere, input: global reference frame
[azimuth,elevation,r] = cart2sph(ESS0(:,1),ESS0(:,2),ESS0(:,3)-zc*m2mm); % The z-coordinate is translated; azimuth and elevation are in rad, r is in mm
PhiTheta = deg2rad(azel2phitheta(rad2deg([azimuth, elevation]'),false)); % First row is phi, second row is theta
Npoints = size(PhiTheta, 2);
SaC.azimuth = azimuth;
SaC.elevation = elevation;

if abs(Ra*m2mm - r) > 1e-10 % check for errors in the translation and in the conversion
    warning("The radius should be constant and equal to Ra")
end

% From Cartesian to spherical components (reference frame of the sphere)
vEs_num = zeros(3,Npoints);
vHs_num = zeros(3,Npoints);
for idx = 1:Npoints
    % Input angles must be in degrees! [az, el, r]
    vEs_num(:,idx) = cart2sphvec(transpose([Fields_SS0.Ex(idx), Fields_SS0.Ey(idx), Fields_SS0.Ez(idx)]),rad2deg(azimuth(idx)),rad2deg(elevation(idx))); % azimuth, elevation, z
    vHs_num(:,idx) = cart2sphvec(transpose([Fields_SS0.Hx(idx), Fields_SS0.Hy(idx), Fields_SS0.Hz(idx)]),rad2deg(azimuth(idx)),rad2deg(elevation(idx)));
end
[vEs_ph, vEs_th, vEs_r, vHs_ph, vHs_th, vHs_r] = azelr2thphr(vEs_num, vHs_num);

%% Integration

% Multipole indices
% Nl = 15; % It should be L
Lmax = 20; % It should be L
% N_idx = 0.5*(Nl+1)*(Nl+2) - 1; % Only positive indexes
N_idx = (Lmax+1)^2 - 1; % Also negative indexes
idx_lm = zeros(N_idx, 2);
idx = 1;
for l = 1:Lmax
    for m = -l:l
        idx_lm(idx, :) = [l, m];
        idx = idx+1;
    end
end
idx_lm_lin = 1:N_idx;
if N_idx ~= size(idx_lm, 1)
    warning("There is an error in the number of indices")
end

SaC.typ = "h2";
[I_Anm, I_BZnm, I_Anm_R, I_BZnm_R] = MultipoleCoefficients2(PhiTheta, vEs_ph, vEs_th, vEs_r, vHs_r, Ra, SaC, N_ph, N_th, N_idx, Npoints, idx_lm);

disp("Maximum relative error from the two ways to compute the coefficients")
disp([max(abs((I_Anm_R - I_Anm)./I_Anm)), max(abs((I_BZnm_R - I_BZnm)./I_BZnm))])
% hh1 = [I_Anm_R,   I_Anm]; hh1(1:20, :)
% hh2 = [I_BZnm_R, I_BZnm]; hh2(1:20, :)

if nnz(isnan(I_Anm)) || nnz(isnan(I_BZnm))
    warning("There is an error in the computation of the multipoles")
end

%% Field reconstruction - method A
% Reconstructed field components in the reference frame of the sphere

SaC.typ = "h2";
[E_sph, H_sph, vEr_Cc, vHr_Cc] = EHfieldME(I_Anm, I_BZnm, idx_lm, Ra, PhiTheta, SaC);
% [E_sph, H_sph, vEr_Cc, vHr_Cc] = EHfieldME(I_Anm_R, I_BZnm_R, idx_lm, Ra, PhiTheta, SaC);

[thTest, phTest, rTest] = CoorTran(C0(:) + [0;Rs;0]);
[E_sph_PL, H_sph_PL, vEr_Cc_PL, vHr_Cc_PL] = EHfieldME(I_Anm, I_BZnm, idx_lm, rTest, [phTest; thTest], SaC);

%% Figures for the reconstructed fields

if abs(SphereIntegral(ones(size(PhiTheta,2),1), PhiTheta, N_ph, Npoints) - 4*pi)/(4*pi) > 1e-3
    warning("There is an error")
end

Rtest = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]*Ra; % , 1000000
N_Rtest = numel(Rtest);

% Flux of the Poynting vector
Spar = zeros(2,N_Rtest,N_idx);
for idx_i = 1:N_Rtest
    for idx_j = 1:N_idx
        n_spec = idx_lm(idx_j,1);
        Spar(1,idx_i,idx_j) = Rtest(idx_i)^2 * 0.5*conj(1j/Z0)*abs( I_Anm(idx_j))^2*n_spec*(n_spec+1)*(-sph_factor(n_spec,k0,Rtest(idx_i),"h2")*(conj(sph_z(n_spec,k0,Rtest(idx_i),"h2"))));
        Spar(2,idx_i,idx_j) = Rtest(idx_i)^2 * 0.5*conj(1j/Z0)*abs(I_BZnm(idx_j))^2*n_spec*(n_spec+1)*(conj(-sph_factor(n_spec,k0,Rtest(idx_i),"h2"))*(-sph_z(n_spec,k0,Rtest(idx_i),"h2")));
    end
end
if nnz(real(Spar) < 0)
    warning("There is an error")
end
Integral_check4_r = transpose(sum(sum(Spar,3),1));

SaC.typ = "h2";
Flux = FluxSMPE(I_Anm, I_BZnm, idx_lm, Rtest(1), Z0, k0, SaC);

Prad_CST = 0.4964235; % Radiated power according to CST (W)
Integral_check4_r_par = zeros(N_Rtest, Lmax);
for index = 1:Lmax
    n_spec = index;
    Integral_check4_r_par(:,index) = transpose(sum(sum(Spar(:,:,idx_lm(:,1)<=n_spec),3),1));
end

Integral_check4_r_par_relError = 100*abs(real(Integral_check4_r_par) - real(Integral_check4_r))./real(Integral_check4_r);
Integral_check4_r_par_relError_CST = 100*abs(real(Integral_check4_r_par) - Prad_CST)./Prad_CST;

RelDiffSteps = 100*abs(diff(Integral_check4_r_par(1,:))./Integral_check4_r_par(1,2:end));

disp("Flux of the Poynting vector on the sphere enclosing the antenna (W)")
disp(Integral_check4_r(1))

%% Multipole coefficients on the scatterer surface from the simulated fields

[azimuth_SSS0,elevation_SSS0,r_SSS0] = cart2sph(ESSS0(:,1)-x0*m2mm,ESSS0(:,2)-y0*m2mm,ESSS0(:,3)-(z0+zc)*m2mm); % The z-coordinate is translated; azimuth and elevation are in rad, r is in mm
PhiTheta_SSS0 = deg2rad(azel2phitheta(rad2deg([azimuth_SSS0, elevation_SSS0]'),false)); % First row is phi, second row is theta
Npoints_SSS0 = size(PhiTheta_SSS0, 2);

if abs(Rs*m2mm - r_SSS0) > 1e-10 % check for errors in the translation and in the conversion
    warning("The radius should be constant and equal to Rs")
end

% max(abs(PhiTheta_SSS0 - PhiTheta), [], 'all')
% min(abs(PhiTheta_SSS0), [], "all")

% From Cartesian to spherical components (reference frame of the sphere)
vEs_num_SSS0 = zeros(3,Npoints_SSS0);
vHs_num_SSS0 = zeros(3,Npoints_SSS0);
for idx = 1:Npoints_SSS0
    % Input angles must be in degrees! [az, el, r]
    vEs_num_SSS0(:,idx) = cart2sphvec(transpose([Fields_SSS0.Ex(idx), Fields_SSS0.Ey(idx), Fields_SSS0.Ez(idx)]),rad2deg(azimuth_SSS0(idx)),rad2deg(elevation_SSS0(idx))); % azimuth, elevation, z
    vHs_num_SSS0(:,idx) = cart2sphvec(transpose([Fields_SSS0.Hx(idx), Fields_SSS0.Hy(idx), Fields_SSS0.Hz(idx)]),rad2deg(azimuth_SSS0(idx)),rad2deg(elevation_SSS0(idx)));
end
[vEs_ph_SSS0, vEs_th_SSS0, vEs_r_SSS0, vHs_ph_SSS0, vHs_th_SSS0, vHs_r_SSS0] = azelr2thphr(vEs_num_SSS0, vHs_num_SSS0);

SaC.typ = "j";
[I_Anm_SSS0, I_BZnm_SSS0, I_Anm_R_SSS0, I_BZnm_R_SSS0] = MultipoleCoefficients2(PhiTheta, vEs_ph_SSS0, vEs_th_SSS0, vEs_r_SSS0, vHs_r_SSS0, Rs, SaC, N_ph, N_th, N_idx, Npoints_SSS0, idx_lm);

disp("Maximum relative error from the two ways to compute the coefficients")
disp([max(abs((I_Anm_R_SSS0 - I_Anm_SSS0)./I_Anm_SSS0)), max(abs((I_BZnm_R_SSS0 - I_BZnm_SSS0)./I_BZnm_SSS0))])
% hh1 = [I_Anm_R,   I_Anm]; hh1(1:20, :)
% hh2 = [I_BZnm_R, I_BZnm]; hh2(1:20, :)

if nnz(isnan(I_Anm_SSS0)) || nnz(isnan(I_BZnm_SSS0))
    warning("There is an error in the computation of the multipoles")
end

%% Transform of the coordinate system - method A

LmaxAT = 15; % From the convergence of the multipoles
idxAT = idx_lm(:,1)<=LmaxAT;

disp("Time elapsed for the translation of the multipole coefficients")
tic
% Check the functions that are used for the radial expansion
SaC.typ = "h2";
[AS, BS, xxx] = CoeffTranslation_v2(x0,y0,z0S,I_Anm(idxAT),I_BZnm(idxAT),k0,idx_lm(idxAT,:),SaC); % z must refer to the reference system of the sphere enclosing the antenna, not to the global reference system
toc

SaC.typ = "j"; % The expansion must be regular in the origin: then, Bessel functions of the first kind must be considered.
[E_sph_S0, H_sph_S0, vEr_S0, vHr_S0] = EHfieldME(AS, BS, idx_lm(idxAT,:), Rs, PhiTheta, SaC);

%% Method based on the Green's function - method B

% PhiTheta
[x,y,z] = sph2cart(azimuth,elevation,Rs);
rO = transpose([x,y,z]); % Observation points in the reference frame of the scatterer
rS = mm2m*transpose([xSides, ySides, zSides]) - transpose(C0); % Source points in the reference frame of the scatterer

Rvec = rS.*ones([size(rS),Npoints]) - permute(rO.*ones([size(rO),size(rS,2)]), [1, 3, 2]);
R = squeeze(sqrt(Rvec(1,:,:).^2 + Rvec(2,:,:).^2 + Rvec(3,:,:).^2));
minRmetB = min(R,[],"all");
maxRmetB = max(R,[],"all");

if minRmetB < mm2m
    warning("The surfaces of the two domains are close")
end

minRmetBapprox = C0(3) - Rs - (Lz/2 + (mean([zmax,zmin]) - zc));
maxRmetBapprox = minRmetBapprox + 2*Rs + Lz;

% Fields on the scatterer sphere
for idx = 1:Npoints

    [E_metBp, H_metBp] = IntMethodB(rO(:,idx), rS, JeR, JmR, Sidx, Sdd, ERS0, SaC);

    E_metB.x(:,idx) = E_metBp.x;
    E_metB.y(:,idx) = E_metBp.y;
    E_metB.z(:,idx) = E_metBp.z;

    % There is an error in H_metBp
    H_metB.x(:,idx) = H_metBp.x;
    H_metB.y(:,idx) = H_metBp.y;
    H_metB.z(:,idx) = H_metBp.z;

end

vEs_MB = zeros(3,Npoints);
vHs_MB = zeros(3,Npoints);
for idx = 1:Npoints
    % Input angles must be in degrees!
    vEs_MB(:,idx) = cart2sphvec(transpose([E_metB.x(idx), E_metB.y(idx), E_metB.z(idx)]),rad2deg(azimuth(idx)),rad2deg(elevation(idx))); % azimuth, elevation, z
    vHs_MB(:,idx) = cart2sphvec(transpose([H_metB.x(idx), H_metB.y(idx), H_metB.z(idx)]),rad2deg(azimuth(idx)),rad2deg(elevation(idx)));
end
[vEs_ph_MB, vEs_th_MB, vEs_r_MB, vHs_ph_MB, vHs_th_MB, vHs_r_MB] = azelr2thphr(vEs_MB, vHs_MB);

% Multipole coefficients in the reference frame of the scatterer
SaC.typ = "j";
[I_Anm_MB, I_BZnm_MB, I_Anm_R_MB, I_BZnm_R_MB] = MultipoleCoefficients2(PhiTheta, vEs_ph_MB, vEs_th_MB, vEs_r_MB, vHs_r_MB, Rs, SaC, N_ph, N_th, N_idx, Npoints, idx_lm);

SaC.typ = "j"; % The expansion must be regular in the origin: then, Bessel functions of the first kind must be considered.
[E_sph_S0_metB, H_sph_S0_metB, vEr_metB, vHr_metB] = EHfieldME(I_Anm_MB, I_BZnm_MB, idx_lm, Rs, PhiTheta, SaC);

SaC.typ = "j";
Flux = FluxSMPE(I_Anm_MB, I_BZnm_MB, idx_lm, Rs, Z0, k0, SaC);

%% Figures

% Compare the fields from all the methods on the scatterer surface

if 0
    run(figures_v2.m)
end

%% Functions

function [vEs_ph, vEs_th, vEs_r, vHs_ph, vHs_th, vHs_r] = azelr2thphr(vEs_num, vHs_num)

% Input:  az, el, r
% Output: ph, th, r

vEs_ph = vEs_num(1, :); vEs_th = -vEs_num(2, :); vEs_r  = vEs_num(3, :);
vHs_ph = vHs_num(1, :); vHs_th = -vHs_num(2, :); vHs_r  = vHs_num(3, :);

end

function [I_Anm, I_BZnm, I_Anm_R, I_BZnm_R] = MultipoleCoefficients2(PhiTheta, vEs_ph, vEs_th, vEs_r, vHs_r, R, SaC, N_ph, N_th, N_idx, Npoints, idx_lm)
% Computation of the multipole coefficients
% They are computed based on
% 1) the tangential components of the electric field
% 2) the radial components of the electric and magnetic fields

k0 = SaC.k0;

% Computation of TSMF
m_th = zeros(N_idx, Npoints);
m_ph = zeros(N_idx, Npoints);
Ylm  = zeros(N_idx, Npoints);
coeff_A = zeros(N_idx, 1);
coeff_B = zeros(N_idx, 1);
coeff_R = zeros(N_idx, 1);

% Computation of scalar harmonics and transverse spherical multipole functions
for idx_lm_spec = 1:N_idx
    % disp([idx_lm(idx_lm_spec, 1), idx_lm(idx_lm_spec, 2)])

    l_spec = idx_lm(idx_lm_spec, 1);
    m_spec = idx_lm(idx_lm_spec, 2);
    [m_th(idx_lm_spec, :), m_ph(idx_lm_spec, :)] = TSMF(PhiTheta(2, :), PhiTheta(1, :), l_spec, m_spec);

    Ylm(idx_lm_spec, :) = sph_harmonC(l_spec, m_spec, PhiTheta(2, :), PhiTheta(1, :));

    % Coefficients for 1)
    coeff_A(idx_lm_spec) = -1./(l_spec*(l_spec+1)*sph_factor(l_spec, k0, R, SaC.typ));
    coeff_B(idx_lm_spec) =  1./(l_spec*(l_spec+1)*sph_z(l_spec, k0, R, SaC.typ));

    % Coefficients for 2)
    coeff_R(idx_lm_spec) = -(k0*R)./(l_spec*(l_spec+1)*sph_z(l_spec, k0, R, SaC.typ));
end

if nnz(isnan(m_th)) || nnz(isnan(m_ph))
    warning("There is an error in the computation of the TSMF")
end

n_th = m_ph;
n_ph =-m_th;

% spy(m_theta==0)
% nnz(m_theta(:,1)==0)

ph_values = PhiTheta(1,1:N_ph); % N or NN, according to the algorithm for the points
% Phi step
d_phi_help = abs(diff(ph_values));
d_phi = mean(d_phi_help);

th_values = PhiTheta(2, 1:N_ph:end);
if numel(th_values) ~= N_th
    warning("There is an error")
end
% Theta step
d_theta_help = abs(diff(th_values)); % It should be the same value
d_theta = mean(d_theta_help);

% ph_value = 0;

I_BZ_thInt_ph= zeros(N_idx, N_ph);
I_BZ_thInt_th= zeros(N_idx, N_ph);
I_A_thInt_ph = zeros(N_idx, N_ph);
I_A_thInt_th = zeros(N_idx, N_ph);
I_AR_thInt = zeros(N_idx, N_ph);
I_BR_thInt = zeros(N_idx, N_ph);
for idx_lm_spec = 1:N_idx

        % integration in theta from 0 to pi
        for idx = 1:N_ph

            ph_value = ph_values(idx);
            idx_angles = idx : N_ph : Npoints; % PhiTheta(1,:) == ph_value; % phi is fixed, theta varies
            if nnz(idx_angles) ~= N_ph/2
                warning("There is an error")
            end

            % zero padding for th = 0 and th = pi
            vs_ph_values = [0, vEs_ph(idx_angles), 0];
            vs_th_values = [0, vEs_th(idx_angles), 0];
            
            TSMF_values_m_ph = [0, m_ph(idx_lm_spec, idx_angles), 0];
            TSMF_values_m_th = [0, m_th(idx_lm_spec, idx_angles), 0];
            TSMF_values_n_ph = [0, n_ph(idx_lm_spec, idx_angles), 0];
            TSMF_values_n_th = [0, n_th(idx_lm_spec, idx_angles), 0];

            sin_values = [0, sin(PhiTheta(2, idx_angles)), 0];

            % Coefficients from the transverse components
            I_A_thInt_ph(idx_lm_spec, idx) = d_theta*trapz(vs_ph_values .* conj(TSMF_values_n_ph) .* sin_values);
            I_A_thInt_th(idx_lm_spec, idx) = d_theta*trapz(vs_th_values .* conj(TSMF_values_n_th) .* sin_values);

            I_BZ_thInt_ph(idx_lm_spec, idx)= d_theta*trapz(vs_ph_values .* conj(TSMF_values_m_ph) .* sin_values);
            I_BZ_thInt_th(idx_lm_spec, idx)= d_theta*trapz(vs_th_values .* conj(TSMF_values_m_th) .* sin_values);

            % Coefficients from the radial components
            vEs_R_values = [0, vEs_r(idx_angles), 0];
            vHs_R_values = [0, vHs_r(idx_angles), 0];

            Ylm_values   = [0, Ylm(idx_lm_spec, idx_angles), 0];

            I_AR_thInt(idx_lm_spec, idx) = d_theta*trapz(vEs_R_values .* conj(Ylm_values) .* sin_values);
            I_BR_thInt(idx_lm_spec, idx) = d_theta*trapz(vHs_R_values .* conj(Ylm_values) .* sin_values);

        end
end
if nnz(isnan(I_BZ_thInt_ph)) || nnz(isnan(I_BZ_thInt_th))
    warning("There is something wrong in the integration along theta")
end

% integration in phi from 0 to 2*pi extension for phi = 2*pi
ph_values_int = [ph_values, 2*pi];

% Coefficients from the transverse components
I_A_thInt_ph = [I_A_thInt_ph, I_A_thInt_ph(:, 1)];
I_A_thInt_th = [I_A_thInt_th, I_A_thInt_th(:, 1)];

I_BZ_thInt_ph = [I_BZ_thInt_ph, I_BZ_thInt_ph(:, 1)];
I_BZ_thInt_th = [I_BZ_thInt_th, I_BZ_thInt_th(:, 1)];

I_A_ph = d_phi*trapz(I_A_thInt_ph, 2);
I_A_th = d_phi*trapz(I_A_thInt_th, 2);

I_BZ_ph = d_phi*trapz(I_BZ_thInt_ph, 2);
I_BZ_th = d_phi*trapz(I_BZ_thInt_th, 2);

% Multipole coefficients
I_BZnm= coeff_B .* (I_BZ_ph + I_BZ_th);
I_Anm = coeff_A .* (I_A_ph + I_A_th);

% Coefficients from the radial components
I_AR = d_phi*trapz([I_AR_thInt, I_AR_thInt(:, 1)], 2);
I_BR = d_phi*trapz([I_BR_thInt, I_BR_thInt(:, 1)], 2);

% Multipole coefficients
I_Anm_R = coeff_R .* I_AR; % I do not know where the minus comes from
I_BZnm_R= coeff_R .* I_BR * (SaC.Z0/(1i));

end

function [I_Anm, I_BZnm, I_Anm_R, I_BZnm_R] = MultipoleCoefficients1(PhiTheta, vEs_ph, vEs_th, vEs_r, vHs_r, R, SaC, N_ph, N_th, N_idx, Npoints, idx_lm)
% Computation of the multipole coefficients
% They are computed based on
% 1) the tangential components of the electric field
% 2) the radial components of the electric and magnetic fields

k0 = SaC.k0;

% Computation of TSMF
m_th = zeros(N_idx, Npoints);
m_ph = zeros(N_idx, Npoints);
Ylm  = zeros(N_idx, Npoints);
coeff_A = zeros(N_idx, 1);
coeff_B = zeros(N_idx, 1);
coeff_R = zeros(N_idx, 1);

% Computation of scalar harmonics and transverse spherical multipole functions
for idx_lm_spec = 1:N_idx
    % disp([idx_lm(idx_lm_spec, 1), idx_lm(idx_lm_spec, 2)])

    l_spec = idx_lm(idx_lm_spec, 1);
    m_spec = idx_lm(idx_lm_spec, 2);
    [m_th(idx_lm_spec, :), m_ph(idx_lm_spec, :)] = TSMF(PhiTheta(2, :), PhiTheta(1, :), l_spec, m_spec);

    Ylm(idx_lm_spec, :) = sph_harmonC(l_spec, m_spec, PhiTheta(2, :), PhiTheta(1, :));

    % Coefficients for 1)
    coeff_A(idx_lm_spec) = -1./(l_spec*(l_spec+1)*sph_factor(l_spec, k0, R, SaC.typ));
    coeff_B(idx_lm_spec) =  1./(l_spec*(l_spec+1)*sph_z(l_spec, k0, R, SaC.typ));

    % Coefficients for 2)
    coeff_R(idx_lm_spec) = -(k0*R)./(l_spec*(l_spec+1)*sph_z(l_spec, k0, R, SaC.typ));
end

if nnz(isnan(m_th)) || nnz(isnan(m_ph))
    warning("There is an error in the computation of the TSMF")
end

n_th = m_ph;
n_ph =-m_th;

% spy(m_theta==0)
% nnz(m_theta(:,1)==0)

ph_values = PhiTheta(1,1:N_ph); % N or NN, according to the algorithm for the points
% Phi step
d_phi_help = abs(diff(ph_values));
d_phi = mean(d_phi_help);

th_values = PhiTheta(2, 1:N_ph:end);
if numel(th_values) ~= N_th
    warning("There is an error")
end
% Theta step
d_theta_help = abs(diff(th_values)); % It should be the same value
d_theta = mean(d_theta_help);

% ph_value = 0;

I_BZ_thInt_ph= zeros(N_idx, N_ph);
I_BZ_thInt_th= zeros(N_idx, N_ph);
I_A_thInt_ph = zeros(N_idx, N_ph);
I_A_thInt_th = zeros(N_idx, N_ph);
I_AR_thInt = zeros(N_idx, N_ph);
I_BR_thInt = zeros(N_idx, N_ph);
for idx_lm_spec = 1:N_idx

        % integration in theta from 0 to pi
        for idx = 1:N_ph

            ph_value = ph_values(idx);
            idx_angles = idx : N_ph : Npoints; % PhiTheta(1,:) == ph_value; % phi is fixed, theta varies
            if nnz(idx_angles) ~= N_ph/2
                warning("There is an error")
            end

            % zero padding for th = 0 and th = pi
            vs_ph_values = [0, vEs_ph(idx_angles), 0];
            vs_th_values = [0, vEs_th(idx_angles), 0];
            
            TSMF_values_m_ph = [0, m_ph(idx_lm_spec, idx_angles), 0];
            TSMF_values_m_th = [0, m_th(idx_lm_spec, idx_angles), 0];
            TSMF_values_n_ph = [0, n_ph(idx_lm_spec, idx_angles), 0];
            TSMF_values_n_th = [0, n_th(idx_lm_spec, idx_angles), 0];

            if issorted(PhiTheta(2, idx_angles))
                th_values_spec = [0, PhiTheta(2, idx_angles), pi];
            else
                th_values_spec = [pi, PhiTheta(2, idx_angles), 0];
                warning("Theta angles are reversed")
            end
            sin_values = sin(th_values_spec);

            % Coefficients from the transverse components
            I_A_thInt_ph(idx_lm_spec, idx) = trapz(th_values_spec, vs_ph_values .* conj(TSMF_values_n_ph) .* sin_values);
            I_A_thInt_th(idx_lm_spec, idx) = trapz(th_values_spec, vs_th_values .* conj(TSMF_values_n_th) .* sin_values);

            I_BZ_thInt_ph(idx_lm_spec, idx)= trapz(th_values_spec, vs_ph_values .* conj(TSMF_values_m_ph) .* sin_values);
            I_BZ_thInt_th(idx_lm_spec, idx)= trapz(th_values_spec, vs_th_values .* conj(TSMF_values_m_th) .* sin_values);

            % Coefficients from the radial components
            vEs_R_values = [0, vEs_r(idx_angles), 0];
            vHs_R_values = [0, vHs_r(idx_angles), 0];

            Ylm_values   = [0, Ylm(idx_lm_spec, idx_angles), 0];

            I_AR_thInt(idx_lm_spec, idx) = trapz(th_values_spec, vEs_R_values .* conj(Ylm_values) .* sin_values);
            I_BR_thInt(idx_lm_spec, idx) = trapz(th_values_spec, vHs_R_values .* conj(Ylm_values) .* sin_values);

        end
end
if nnz(isnan(I_BZ_thInt_ph)) || nnz(isnan(I_BZ_thInt_th))
    warning("There is something wrong in the integration along theta")
end

% integration in phi from 0 to 2*pi extension for phi = 2*pi
ph_values_int = [ph_values, 2*pi];

% Coefficients from the transverse components
I_A_thInt_ph = [I_A_thInt_ph, I_A_thInt_ph(:, 1)];
I_A_thInt_th = [I_A_thInt_th, I_A_thInt_th(:, 1)];

I_BZ_thInt_ph = [I_BZ_thInt_ph, I_BZ_thInt_ph(:, 1)];
I_BZ_thInt_th = [I_BZ_thInt_th, I_BZ_thInt_th(:, 1)];

I_A_ph = trapz(ph_values_int, I_A_thInt_ph, 2);
I_A_th = trapz(ph_values_int, I_A_thInt_th, 2);

I_BZ_ph = trapz(ph_values_int, I_BZ_thInt_ph, 2);
I_BZ_th = trapz(ph_values_int, I_BZ_thInt_th, 2);

% Multipole coefficients
I_BZnm= coeff_B .* (I_BZ_ph + I_BZ_th);
I_Anm = coeff_A .* (I_A_ph + I_A_th);

% Coefficients from the radial components
I_AR = trapz(ph_values_int, [I_AR_thInt, I_AR_thInt(:, 1)], 2);
I_BR = trapz(ph_values_int, [I_BR_thInt, I_BR_thInt(:, 1)], 2);

% Multipole coefficients
I_Anm_R = coeff_R .* I_AR; % I do not know where the minus comes from
I_BZnm_R= coeff_R .* I_BR * (SaC.Z0/(1i));

end

function [Je, Jm, JeR, JmR, xSides, ySides, zSides] = EquivalentCurrents(ERS0, Fields_RS0, Sidx)

% Je: equivalent electric current as cell array
% Jm: equivalent magnetic current as cell array
% JeR: equivalent electric current as struct
% JmR: equivalent magnetic current as struct
% xSides, ySides, zSides: x, y, and z coordinates of the points

% The phasors must be used (they are complex numbers)
Je{1}.x = zeros(nnz(Sidx{1}),1);
Je{1}.y = +Fields_RS0.Hz(Sidx{1});
Je{1}.z = -Fields_RS0.Hy(Sidx{1});

Je{2}.x = zeros(nnz(Sidx{2}),1);
Je{2}.y = -Fields_RS0.Hz(Sidx{2});
Je{2}.z = +Fields_RS0.Hy(Sidx{2});

Je{3}.x = -Fields_RS0.Hz(Sidx{3});
Je{3}.y = zeros(nnz(Sidx{3}),1);
Je{3}.z = +Fields_RS0.Hx(Sidx{3});

Je{4}.x = +Fields_RS0.Hz(Sidx{4});
Je{4}.y = zeros(nnz(Sidx{4}),1);
Je{4}.z = -Fields_RS0.Hx(Sidx{4});

Je{5}.x = +Fields_RS0.Hy(Sidx{5});
Je{5}.y = -Fields_RS0.Hx(Sidx{5});
Je{5}.z = zeros(nnz(Sidx{5}),1);

Je{6}.x = -Fields_RS0.Hy(Sidx{6});
Je{6}.y = +Fields_RS0.Hx(Sidx{6});
Je{6}.z = zeros(nnz(Sidx{6}),1);

Jm{1}.x = zeros(nnz(Sidx{1}),1);
Jm{1}.y = -Fields_RS0.Ez(Sidx{1});
Jm{1}.z = +Fields_RS0.Ey(Sidx{1});

Jm{2}.x = zeros(nnz(Sidx{2}),1);
Jm{2}.y = +Fields_RS0.Ez(Sidx{2});
Jm{2}.z = -Fields_RS0.Ey(Sidx{2});

Jm{3}.x = +Fields_RS0.Ez(Sidx{3});
Jm{3}.y = zeros(nnz(Sidx{3}),1);
Jm{3}.z = -Fields_RS0.Ex(Sidx{3});

Jm{4}.x = -Fields_RS0.Ez(Sidx{4});
Jm{4}.y = zeros(nnz(Sidx{4}),1);
Jm{4}.z = +Fields_RS0.Ex(Sidx{4});

Jm{5}.x = -Fields_RS0.Ey(Sidx{5});
Jm{5}.y = +Fields_RS0.Ex(Sidx{5});
Jm{5}.z = zeros(nnz(Sidx{5}),1);

Jm{6}.x = +Fields_RS0.Ey(Sidx{6});
Jm{6}.y = -Fields_RS0.Ex(Sidx{6});
Jm{6}.z = zeros(nnz(Sidx{6}),1);

xSides = ERS0(Sidx{1},1);
ySides = ERS0(Sidx{1},2);
zSides = ERS0(Sidx{1},3);
JeR.x = Je{1}.x;
JeR.y = Je{1}.y;
JeR.z = Je{1}.z;
JmR.x = Jm{1}.x;
JmR.y = Jm{1}.y;
JmR.z = Jm{1}.z;
for idx = 2:6
    xSides = [xSides; ERS0(Sidx{idx},1)];
    ySides = [ySides; ERS0(Sidx{idx},2)];
    zSides = [zSides; ERS0(Sidx{idx},3)];

    JeR.x = [JeR.x; Je{idx}.x];
    JeR.y = [JeR.y; Je{idx}.y];
    JeR.z = [JeR.z; Je{idx}.z];

    JmR.x = [JmR.x; Jm{idx}.x];
    JmR.y = [JmR.y; Jm{idx}.y];
    JmR.z = [JmR.z; Jm{idx}.z];
end

end

function [E, H] = IntMethodB(rO, rS, JeR, JmR, Sidx, Sdd, ERS0, SaC)

k0 = SaC.k0;
Z0 = SaC.Z0;
Y0 = SaC.Y0;

h1 = DGF(k0, rO, rS);
h2 = GSGFcross(k0, rO, rS);

comp_x_E_Je = (squeeze(h1(1,1,:)) .* JeR.x + squeeze(h1(1,2,:)) .* JeR.y + squeeze(h1(1,3,:)) .* JeR.z);
comp_y_E_Je = (squeeze(h1(2,1,:)) .* JeR.x + squeeze(h1(2,2,:)) .* JeR.y + squeeze(h1(2,3,:)) .* JeR.z);
comp_z_E_Je = (squeeze(h1(3,1,:)) .* JeR.x + squeeze(h1(3,2,:)) .* JeR.y + squeeze(h1(3,3,:)) .* JeR.z);

comp_x_E_Jm = (squeeze(h2(1,1,:)) .* JmR.x + squeeze(h2(1,2,:)) .* JmR.y + squeeze(h2(1,3,:)) .* JmR.z);
comp_y_E_Jm = (squeeze(h2(2,1,:)) .* JmR.x + squeeze(h2(2,2,:)) .* JmR.y + squeeze(h2(2,3,:)) .* JmR.z);
comp_z_E_Jm = (squeeze(h2(3,1,:)) .* JmR.x + squeeze(h2(3,2,:)) .* JmR.y + squeeze(h2(3,3,:)) .* JmR.z);

comp_x_H_Je = (squeeze(h2(1,1,:)) .* JeR.x + squeeze(h2(1,2,:)) .* JeR.y + squeeze(h2(1,3,:)) .* JeR.z);
comp_y_H_Je = (squeeze(h2(2,1,:)) .* JeR.x + squeeze(h2(2,2,:)) .* JeR.y + squeeze(h2(2,3,:)) .* JeR.z);
comp_z_H_Je = (squeeze(h2(3,1,:)) .* JeR.x + squeeze(h2(3,2,:)) .* JeR.y + squeeze(h2(3,3,:)) .* JeR.z);

comp_x_H_Jm = (squeeze(h1(1,1,:)) .* JmR.x + squeeze(h1(1,2,:)) .* JmR.y + squeeze(h1(1,3,:)) .* JmR.z);
comp_y_H_Jm = (squeeze(h1(2,1,:)) .* JmR.x + squeeze(h1(2,2,:)) .* JmR.y + squeeze(h1(2,3,:)) .* JmR.z);
comp_z_H_Jm = (squeeze(h1(3,1,:)) .* JmR.x + squeeze(h1(3,2,:)) .* JmR.y + squeeze(h1(3,3,:)) .* JmR.z);

I_comp_x_E_Je = 0; I_comp_y_E_Je = 0; I_comp_z_E_Je = 0;
I_comp_x_E_Jm = 0; I_comp_y_E_Jm = 0; I_comp_z_E_Jm = 0;
I_comp_x_H_Je = 0; I_comp_y_H_Je = 0; I_comp_z_H_Je = 0;
I_comp_x_H_Jm = 0; I_comp_y_H_Jm = 0; I_comp_z_H_Jm = 0;
for idx = 1:6 % Integrals over the six faces of the enclosing parallelepiped

    weights = ones(nnz(Sidx{idx}),1); % 1/4 on the corners, 1/2 on the sides, 1 in the middle
    if idx == 1 || idx == 2
        weights(ERS0(Sidx{idx},2) == SaC.ymaxSim) = 0.5*weights(ERS0(Sidx{idx},2) == SaC.ymaxSim);
        weights(ERS0(Sidx{idx},2) == SaC.yminSim) = 0.5*weights(ERS0(Sidx{idx},2) == SaC.yminSim);
        weights(ERS0(Sidx{idx},3) == SaC.zmaxSim) = 0.5*weights(ERS0(Sidx{idx},3) == SaC.zmaxSim);
        weights(ERS0(Sidx{idx},3) == SaC.zminSim) = 0.5*weights(ERS0(Sidx{idx},3) == SaC.zminSim);
    elseif idx == 3 || idx == 4
        weights(ERS0(Sidx{idx},1) == SaC.xmaxSim) = 0.5*weights(ERS0(Sidx{idx},1) == SaC.xmaxSim);
        weights(ERS0(Sidx{idx},1) == SaC.xminSim) = 0.5*weights(ERS0(Sidx{idx},1) == SaC.xminSim);
        weights(ERS0(Sidx{idx},3) == SaC.zmaxSim) = 0.5*weights(ERS0(Sidx{idx},3) == SaC.zmaxSim);
        weights(ERS0(Sidx{idx},3) == SaC.zminSim) = 0.5*weights(ERS0(Sidx{idx},3) == SaC.zminSim);
    elseif idx == 5 || idx == 6
        weights(ERS0(Sidx{idx},1) == SaC.xmaxSim) = 0.5*weights(ERS0(Sidx{idx},1) == SaC.xmaxSim);
        weights(ERS0(Sidx{idx},1) == SaC.xminSim) = 0.5*weights(ERS0(Sidx{idx},1) == SaC.xminSim);
        weights(ERS0(Sidx{idx},2) == SaC.ymaxSim) = 0.5*weights(ERS0(Sidx{idx},2) == SaC.ymaxSim);
        weights(ERS0(Sidx{idx},2) == SaC.yminSim) = 0.5*weights(ERS0(Sidx{idx},2) == SaC.yminSim);
    end

    if nnz(weights==0.25) ~= 4
        warning("There is an error")
    end

    I_comp_x_E_Je = I_comp_x_E_Je + Sdd(idx)*sum(weights.*comp_x_E_Je(Sidx{idx}), "all");
    I_comp_y_E_Je = I_comp_y_E_Je + Sdd(idx)*sum(weights.*comp_y_E_Je(Sidx{idx}), "all");
    I_comp_z_E_Je = I_comp_z_E_Je + Sdd(idx)*sum(weights.*comp_z_E_Je(Sidx{idx}), "all");

    I_comp_x_E_Jm = I_comp_x_E_Jm + Sdd(idx)*sum(weights.*comp_x_E_Jm(Sidx{idx}), "all");
    I_comp_y_E_Jm = I_comp_y_E_Jm + Sdd(idx)*sum(weights.*comp_y_E_Jm(Sidx{idx}), "all");
    I_comp_z_E_Jm = I_comp_z_E_Jm + Sdd(idx)*sum(weights.*comp_z_E_Jm(Sidx{idx}), "all");

    I_comp_x_H_Je = I_comp_x_H_Je + Sdd(idx)*sum(weights.*comp_x_H_Je(Sidx{idx}), "all");
    I_comp_y_H_Je = I_comp_y_H_Je + Sdd(idx)*sum(weights.*comp_y_H_Je(Sidx{idx}), "all");
    I_comp_z_H_Je = I_comp_z_H_Je + Sdd(idx)*sum(weights.*comp_z_H_Je(Sidx{idx}), "all");

    I_comp_x_H_Jm = I_comp_x_H_Jm + Sdd(idx)*sum(weights.*comp_x_H_Jm(Sidx{idx}), "all");
    I_comp_y_H_Jm = I_comp_y_H_Jm + Sdd(idx)*sum(weights.*comp_y_H_Jm(Sidx{idx}), "all");
    I_comp_z_H_Jm = I_comp_z_H_Jm + Sdd(idx)*sum(weights.*comp_z_H_Jm(Sidx{idx}), "all");
end

E.x = -1j*k0*Z0* I_comp_x_E_Je - I_comp_x_E_Jm;
E.y = -1j*k0*Z0* I_comp_y_E_Je - I_comp_y_E_Jm;
E.z = -1j*k0*Z0* I_comp_z_E_Je - I_comp_z_E_Jm;

H.x = -1j*k0*Y0* I_comp_x_H_Jm + I_comp_x_H_Je;
H.y = -1j*k0*Y0* I_comp_y_H_Jm + I_comp_y_H_Je;
H.z = -1j*k0*Y0* I_comp_z_H_Jm + I_comp_z_H_Je;

end

function output = DGF(k, rO, rS)
% k: wavenumber
% rO: observation point in Cartesian coordinates; dimensions: 3xNpoints
% rS: source point in Cartesian coordinates; dimensions: 3xNpoints

Rvec = rO - rS;
R = vecnorm(rO-rS,2);
Rver = Rvec./R;
Rverx = Rver(1,:);
Rvery = Rver(2,:);
Rverz = Rver(3,:);

GF = SGF(k, R);

DGFmagRR = (3.*(k*R).^(-2) + 3j*(k*R).^(-1) - 1) .* GF;
DGFmagI  = (1 - 1j*(k*R).^(-1) - (k*R).^(-2)) .* GF;

tensor = zeros(3,3,numel(R));
tensor(1,1,:) = Rverx.^2;
tensor(1,2,:) = Rverx.*Rvery;
tensor(1,3,:) = Rverx.*Rverz;
tensor(2,1,:) = Rvery.*Rverx;
tensor(2,2,:) = Rvery.^2;
tensor(2,3,:) = Rvery.*Rverz;
tensor(3,1,:) = Rverz.*Rverx;
tensor(3,2,:) = Rverz.*Rvery;
tensor(3,3,:) = Rverz.^2;

aa = permute(repmat(transpose(DGFmagRR),1,3,3), [2, 3, 1]);
bb = permute(repmat(transpose( DGFmagI),1,3,3), [2, 3, 1]);

output1 = aa .* tensor; % dim: 3x3xNpoints
output2 = bb .* [1, 0, 0; 0, 1, 0; 0, 0, 1];

output = output1 + output2;

end

function output = GSGFcross(k, rO, rS)
% k: wavenumber
% rO: observation point in Cartesian coordinates; dimensions: 3xNpoints
% rS: source point in Cartesian coordinates; dimensions: 3xNpoints

Rvec = rO - rS;
R = vecnorm(rO-rS,2);
Rver = Rvec./R;

GF = SGF(k, R);
GSGFmag = (-1j*k - 1./R) .* GF;

GSGFx = GSGFmag .* Rver(1,:); % dim: 1xNpoints
GSGFy = GSGFmag .* Rver(2,:);
GSGFz = GSGFmag .* Rver(3,:);

output = zeros(3,3,numel(R)); % dim: 3x3xNpoints

% output(1,1,:) = 0;
output(1,2,:) = -GSGFz;
output(1,3,:) = GSGFy;
output(2,1,:) = GSGFz;
% output(2,2,:) = 0;
output(2,3,:) = -GSGFx;
output(3,1,:) = -GSGFy;
output(3,2,:) = GSGFx;
% output(3,2,:) = 0;

end

function output = SGF(k, R)
output = exp(-1j*k*R)./(4*pi*R);
end

function Integral = SphereIntegral(integrand, PhiTheta, N_ph, Npoints)

ph_values = PhiTheta(1,1:N_ph);
Integral_theta = zeros(1, N_ph);

for idx = 1:N_ph
            % ph_value = ph_values(idx);
            idx_angles = idx : N_ph : Npoints;
            if nnz(idx_angles) ~= N_ph/2
                warning("There is an error")
            end

            if issorted(PhiTheta(2, idx_angles))
                th_values = [0, PhiTheta(2, idx_angles), pi];
            else
                th_values = [pi, PhiTheta(2, idx_angles), 0];
                warning("Theta angles are reversed")
            end

            sin_values= sin(th_values);
            integrand_values = [0, transpose(integrand(idx_angles)), 0];

%             Integral_theta(idx) = trapz(sin_values, integrand_values .* sin_values );
            Integral_theta(idx) = trapz(th_values, integrand_values .* sin_values );
end

ph_values_int = [ph_values, 2*pi];

Integral = trapz(ph_values_int, [Integral_theta, Integral_theta(end)]);

end

function [CSTfields] = FieldsFromFWS(Efield, Hfield, PaS)

Ex = Efield(:,4) + 1j*Efield(:,5);
Ey = Efield(:,6) + 1j*Efield(:,7);
Ez = Efield(:,8) + 1j*Efield(:,9);

Hx = Hfield(:,4) + 1j*Hfield(:,5);
Hy = Hfield(:,6) + 1j*Hfield(:,7);
Hz = Hfield(:,8) + 1j*Hfield(:,9);

CSTfields.Ex = Ex;
CSTfields.Ey = Ey;
CSTfields.Ez = Ez;

CSTfields.Hx = Hx;
CSTfields.Hy = Hy;
CSTfields.Hz = Hz;

end

function CSTfieldsTimeEvaluation = TimeEvaluation(CSTfields, f0, t)

CSTfieldsTimeEvaluation.Ex = real(CSTfields.Ex * exp(1i*2*pi*f0*t));
CSTfieldsTimeEvaluation.Ey = real(CSTfields.Ey * exp(1i*2*pi*f0*t));
CSTfieldsTimeEvaluation.Ez = real(CSTfields.Ez * exp(1i*2*pi*f0*t));

CSTfieldsTimeEvaluation.Hx = real(CSTfields.Hx * exp(1i*2*pi*f0*t));
CSTfieldsTimeEvaluation.Hy = real(CSTfields.Hy * exp(1i*2*pi*f0*t));
CSTfieldsTimeEvaluation.Hz = real(CSTfields.Hz * exp(1i*2*pi*f0*t));

end

function [E_sph, H_sph, vErect, vHrect] = EHfieldME(A, BZ, idx_lm, Reval, PhiTheta, PaS)
% Field reconstrution from the multipole coefficients
% A and BZ - coefficients of the multipole expansion (V/m)
% idx_lm - multipole coefficients
% Reval - radius where the multipole expansion is evaluated (m)
% PhiTheta - phi and theta angles (rad)

k0 = PaS.k0;
Z0 = PaS.Z0;
azimuth = PaS.azimuth;
elevation = PaS.elevation;

Npoints = size(PhiTheta, 2);
N_idx = size(idx_lm, 1);

M_r  = zeros(N_idx, Npoints);
M_th = zeros(N_idx, Npoints);
M_ph = zeros(N_idx, Npoints);
N_r  = zeros(N_idx, Npoints);
N_th = zeros(N_idx, Npoints);
N_ph = zeros(N_idx, Npoints);

for idx_lm_spec = 1:N_idx

    l_spec = idx_lm(idx_lm_spec, 1);
    m_spec = idx_lm(idx_lm_spec, 2);

    [M_r(idx_lm_spec, :), M_th(idx_lm_spec, :), M_ph(idx_lm_spec, :), N_r(idx_lm_spec, :), N_th(idx_lm_spec, :), N_ph(idx_lm_spec, :)] = VSMF(Reval, PhiTheta(2, :), PhiTheta(1, :), l_spec, m_spec, k0, PaS.typ);
end

E_r  = sum(A .* N_r  + BZ .* M_r );
E_th = sum(A .* N_th + BZ .* M_th);
E_ph = sum(A .* N_ph + BZ .* M_ph);

H_r  = (1j/Z0)*sum(A .* M_r  + BZ .* N_r );
H_th = (1j/Z0)*sum(A .* M_th + BZ .* N_th);
H_ph = (1j/Z0)*sum(A .* M_ph + BZ .* N_ph);

E_sph.E_r  = E_r;
E_sph.E_th = E_th;
E_sph.E_ph = E_ph;

H_sph.H_r  = H_r;
H_sph.H_th = H_th;
H_sph.H_ph = H_ph;

vEsph = [E_ph; -E_th; E_r]; % The minus sign derives from the different convention for the unit vectors
vHsph = [H_ph; -H_th; H_r];

vErect = zeros(3, Npoints);
vHrect = zeros(3, Npoints);
for idxPoints = 1:Npoints
    % Input angles must be in degrees!
    vErect(:,idxPoints) = sph2cartvec(vEsph(:,idxPoints),rad2deg(azimuth(idxPoints)),rad2deg(elevation(idxPoints)));
    vHrect(:,idxPoints) = sph2cartvec(vHsph(:,idxPoints),rad2deg(azimuth(idxPoints)),rad2deg(elevation(idxPoints)));
end

E_Er_S0 = rms(transpose([real(vErect(1,:)); real(vErect(2,:)); real(vErect(3,:))]),"all")^2/Z0* 4*pi*Reval^2;
E_Es_S0 = rms(transpose([real(vEsph(1,:)); real(vEsph(2,:)); real(vEsph(3,:))]),"all")^2/Z0* 4*pi*Reval^2;
if abs(E_Er_S0 - E_Es_S0) > 100*eps
    error("The power-related quantity must be the same in Cartesian and spherical coordinates")
end

if nnz(isnan(E_sph.E_th)) || nnz(isnan(E_sph.E_ph))
    error("There is an error in the reconstruction")
end

end

function PointList(Az, El, Ra, Center, m2mm, N, filename)

[xA,yA,zA] = sph2cart(Az, El, Ra*m2mm); TH = 1e-6;
xAG = xA + Center(1)*m2mm; XAG = reshape(xAG, [],1);
yAG = yA + Center(2)*m2mm; YAG = reshape(yAG, [],1);
zAG = zA + Center(3)*m2mm; ZAG = reshape(zAG, [],1);

XAG_zeros = XAG;
YAG_zeros = YAG;
ZAG_zeros = ZAG;
XAG_zeros(abs(XAG) < TH) = 0;
YAG_zeros(abs(YAG) < TH) = 0;
ZAG_zeros(abs(ZAG) < TH) = 0;

% AntennaSphereCoordinates = [XAG, YAG, ZAG]; % N*(2*N - 1)
AntennaSphereCoordinates = unique([XAG_zeros, YAG_zeros, ZAG_zeros], 'rows', 'stable'); % N*(2*N - 1) - (2*N-2)*2-(N-2)
% if size(AntennaSphereCoordinates,1) ~= N*(2*N - 1) - (2*N-2)*2-(N-2)
%     warning("The threshold is too large")
% end

% figure
% scatter3(XAG,YAG,ZAG)

figure
scatter3(AntennaSphereCoordinates(:,1), AntennaSphereCoordinates(:,2), AntennaSphereCoordinates(:,3))
title("Sampling points")

figure
scatter(AntennaSphereCoordinates(:,1), AntennaSphereCoordinates(:,2))
axis equal
title("Sampling points")

fid = fopen(strcat("FieldData\", filename), 'w+');
fprintf(fid, '%8.3f %8.3f %8.3f\n', transpose(AntennaSphereCoordinates));
fclose(fid);

end

function [N_ph, N_th] = PointList2(Ra, Center, m2mm, L, QH, filename)

[w_gl,x_gl,X_gl] = QpS2(L,QH); % sum(X_gl.^2) % check radius

N_ph = 2*sqrt(size(x_gl,2)/2); % Number of points for phi
N_th = N_ph/2;

X_gl_denorm = X_gl*Ra*m2mm;

X_gl_denorm(1,:) = X_gl_denorm(1,:) + Center(1)*m2mm;
X_gl_denorm(2,:) = X_gl_denorm(2,:) + Center(2)*m2mm;
X_gl_denorm(3,:) = X_gl_denorm(3,:) + Center(3)*m2mm;

figure
scatter3(X_gl_denorm(1,:), X_gl_denorm(2,:), X_gl_denorm(3,:))
title("Sampling points")

fid = fopen(strcat("FieldData\", filename), 'w+');
fprintf(fid, '%8.4f %8.4f %8.4f\n', X_gl_denorm);
fclose(fid);

end

function Flux = FluxSMPE(I_Anm, I_BZnm, idx_lm, Rtest, Z0, k0, SaC)

Lmax = max(idx_lm(:,1));
N_idx = (Lmax+1)^2 - 1;

typ = SaC.typ;

Spar = zeros(2,N_idx);
for idx_j = 1:N_idx
    n_spec = idx_lm(idx_j,1);
    Spar(1,idx_j) = Rtest^2 * 0.5*conj(1j/Z0)*abs( I_Anm(idx_j))^2*n_spec*(n_spec+1)*(-sph_factor(n_spec,k0,Rtest,typ)*(conj(sph_z(n_spec,k0,Rtest,typ))));
    Spar(2,idx_j) = Rtest^2 * 0.5*conj(1j/Z0)*abs(I_BZnm(idx_j))^2*n_spec*(n_spec+1)*(conj(-sph_factor(n_spec,k0,Rtest,typ))*(-sph_z(n_spec,k0,Rtest,typ)));
end

if nnz(real(Spar) < 0)
    warning("There is an error")
end
Flux = transpose(sum(sum(Spar,2),1));

if typ == "j"
    if real(Flux) ~=0
        error("Without sources the flux of S should be purely imaginary")
    else
        disp("Without sources the flux of S is purely imaginary")
    end
end

end

function [th, ph, r] = CoorTran(Pos)
    x = Pos(1,:);
    y = Pos(2,:);
    z = Pos(3,:);
    [azimuth,elevation,r] = cart2sph(x,y,z); % The z-coordinate is translated; azimuth and elevation are in rad, r is in mm
    PhiTheta = deg2rad(azel2phitheta(rad2deg([azimuth, elevation]'),false)); % First row is phi, second row is theta
    ph = PhiTheta(1, 1);
    th = PhiTheta(2, 1);
end
