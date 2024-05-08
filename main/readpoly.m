function poly = readpoly(polyfile)
%=================================================================
% Read polygon coords from text file into cell array
%                                                                  
% INPUT:                                                           
%   polyfile: path to poly file
% OUTPUT:                                                          
%   poly:  cell array of polygons (x y coords), or a single matrix if file
%          only included one poly
%   
% polyfile should be a comma-separated two column text file of x y
% coordinates, with a #-line between polygons.
%
% Andrew Watson @ Leeds, 21/07/2021
%                                                                  
%=================================================================

% open file
tmp = readcell(polyfile);

% find #'s
hash_ind = find(strcmp(tmp,'#'));

% create indexing vectors
ind1 = [1; hash_ind+1];
ind2 = [hash_ind-1; length(tmp)];

% group cells
poly = cell(length(hash_ind)+1,1);
for ii = 1:length(ind1)
    poly{ii} = cell2mat(tmp(ind1(ii):ind2(ii),:));
end

% convert to matrix if only 1 poly
if length(poly) == 1
    poly = cell2mat(poly);
end

end