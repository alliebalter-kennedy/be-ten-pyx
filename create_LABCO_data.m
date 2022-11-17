function [loc he be] = create_LABCO_data(flag)

% This script creates usable data from the LABCO_data.txt file. There is an
% optional flag 'ignore_outliers', which will return only data that are not
% outliers. otherwise, all data are returned. 

if nargin < 1
    flag = '';
end

addpath('data/')

%% location data

l = readtable('LABCO_loc_data.txt');

loc.lat = table2array(l(:, 'lat'));
loc.lon = table2array(l(:, 'lon'));
loc.elv = table2array(l(:, 'elv'));

%% load data

d = readtable('LABCO_data.txt');

%% bring all data into beryllium and helium structures

% beryllium

be.id = table2cell(d(d.nuclide == 10, 'sample_id'));
be.depth_cm = table2array(d(d.nuclide == 10, 'depth'));
be.aliquot = table2cell(d(d.nuclide == 10, 'aliquot'));
be.prep_lab = table2cell(d(d.nuclide == 10, 'prep_lab'));
be.analysis_lab = table2cell(d(d.nuclide == 10, 'analysis_lab'));
be.outlier = table2array(d(d.nuclide == 10, 'outlier'));

be.N10 = table2array(d(d.nuclide == 10, 'N'));
be.dN10 = table2array(d(d.nuclide == 10, 'dN')); 

% helium

he.id = table2cell(d(d.nuclide == 3, 'sample_id'));
he.depth_cm = table2array(d(d.nuclide == 3, 'depth'));
he.aliquot = table2cell(d(d.nuclide == 3, 'aliquot'));
he.prep_lab = table2cell(d(d.nuclide == 3, 'prep_lab'));
he.analysis_lab = table2cell(d(d.nuclide == 3, 'analysis_lab'));
he.outlier = table2array(d(d.nuclide == 3, 'outlier'));

he.N3 = table2array(d(d.nuclide == 3, 'N'));
he.dN3 = table2array(d(d.nuclide == 3, 'dN')); 
he.N3_normalized = he.N3 .* (5.02E+09 ./ table2array(d(d.nuclide == 3, 'standard_N')));
he.N4 = table2array(d(d.nuclide == 3, 'N4'));
he.dN4 = table2array(d(d.nuclide == 3, 'dN4'));

if strcmp(flag, 'ignore_outliers') == 1
    be.id = be.id(be.outlier == 0, :);
    be.depth_cm = be.depth_cm(be.outlier == 0, :);
    be.aliquot = be.aliquot(be.outlier == 0, :);
    be.prep_lab = be.prep_lab(be.outlier == 0, :);
    be.analysis_lab = be.analysis_lab(be.outlier == 0, :);
    
    be.N10 = be.N10(be.outlier == 0, :);
    be.dN10 = be.dN10(be.outlier == 0, :);
    
    % helium
    
    he.id = he.id(he.outlier == 0, :);
    he.depth_cm = he.depth_cm(he.outlier == 0, :);
    he.aliquot = he.aliquot(he.outlier == 0, :);
    he.prep_lab = he.prep_lab(he.outlier == 0, :);
    he.analysis_lab = he.analysis_lab(he.outlier == 0, :);
    
    he.N3 = he.N3(he.outlier == 0, :);
    he.dN3 = he.dN3(he.outlier == 0, :);
    he.N3_normalized = he.N3_normalized(he.outlier == 0, :);
    he.N4 = he.N4(he.outlier == 0, :);
    he.dN4 = he.dN4(he.outlier == 0, :);


end
end