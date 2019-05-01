
clear;clc;

brain_model = niftiread('brain_model.nii');
volumeViewer(brain_model);

brain_model_mask = brain_model>0;

% solido
brain_fluid_mask = imgaussfilt3(uint8(imdilate(brain_model_mask, ones(9,9,9))), 10)>0;
cranio_mask = imgaussfilt3(uint8(imdilate(brain_fluid_mask,ones(9,9,9))), 10)>0;
pele_mask = imdilate(cranio_mask,ones(5,5,5))>0;

% so borda
pele_mask = pele_mask & ~cranio_mask;
cranio_mask = cranio_mask & ~brain_fluid_mask;
brain_fluid_mask = brain_fluid_mask & ~brain_model_mask;

% adiciona alterações ao modelo
brain_model(brain_fluid_mask)=5;

brain_model(cranio_mask)=255;

brain_model(brain_model>0 & brain_model<21)=21;
brain_model(pele_mask)=15;



volumeViewer(brain_model);
