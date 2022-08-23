function actin_polarity_perform_particle_rec(BasePath,ParticlePath,ParticleName)
%%%%%%%%%
%%%%%%%%%
%
% This function performs particle reconstruction of actin segments.
%
% INPUT
% BasePath --- base path of tiltseries
%              with alignment [BasePath '/alg/alg.mat']
%              and CTF corrected projections [BasePath '/proj/ctfcor']
% Adjust the thickness of tomogram
%
% OUTPUT
% particle_*.em --- particle reconstruction of a actin filament

% Particle reconstruction parameter
Imdim = 4096;%pixel
ot = 1024;%pixel
DimRec = [4096 4096 4096];% /pixel
DimParticle = [160 160 160];% /box size pixel
binning = 0;
tomothickness = 384;

% Display information
disp('Perform particle reconstruction of actin segments from ctfcorrected projections (actin_polarity_cor_final_160.mat)...');
disp(BasePath);

% Load alignment
load([BasePath '/alg/alg.mat'],'algXY');

% Load particle coordinates
load([BasePath '/cor/actin_polarity_cor_final_160.mat'],'plist');

% Align tiltseries
AlignTiltseries(algXY,[BasePath '/proj/ctfcor'],'ctfcor_',[BasePath '/proj/prep'],'prep_');

% Prepare tiltseries
PrepareTiltseries([BasePath '/proj/prep'],'prep_',[BasePath '/proj/prep'],'prep_',size(algXY.Tiltangles,2),algXY.Tiltangles,[],[],ot,Imdim);

% Modify plist
plist(:,3) = plist(:,3) + tomothickness;

% Reconstruct particles
WBPParticlesRec([BasePath '/proj/prep'],'prep_',size(algXY.Tiltangles,2),algXY.Tiltangles,ParticlePath,ParticleName,plist.*4,DimParticle,DimRec,binning,Imdim);

% Delete prepared projections
cd([BasePath '/proj/prep']);
delete *.*;

