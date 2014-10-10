function vpData = sc_read_vpdata(fileName)

% SC_READ_VPDATA: read the data from vanishing point detection algorithm
% Input:
%   - fileName: the txt file containing the pre-computed vanishing point
%   detection code
% Output:
%   - vpData
%   The data structure of vpData
%       - vpData.numVP: number of detected vanishing points
%       - vpData.vp{i}.pos: the vanishing point position in the homogenous coordiante
%       - vpData.vp{i}.score: the score of the vanishing point
%       - vpData.vp{i}.numLines: number of lines supporting the vanishing point
%       - vpData.vp{i}.lines{j}.p1: (x1, y1): starting position
%       - vpData.vp{i}.lines{j}.p2: (x2, y2): ending position
%       - vpData.vp{i}.lines{j}.length: length of the line segment

vpData = [];

% Read data
fid = fopen(fileName);

%% Parse VP positions
temp = fscanf(fid, '%s ', [1 5]);
numVP = 0;
readVPFlag = 1;
VP = [];
while(readVPFlag)
    numVP = numVP + 1;
    vpCurr = fscanf(fid, '%g %g %g %g %g', [5 1]);
    if(~isempty(vpCurr))
        VP(:,numVP) = vpCurr;
    else
        temp = fscanf(fid, '%s ', [1 6]);
        readVPFlag = 0;
    end
end
VP = VP';

vpData.numVP = size(VP, 1);

% Save VP position data
for i = 1: vpData.numVP
    vpData.vp{i}.pos = VP(i, 1:3);
    vpData.vp{i}.score = VP(i, 4);
    vpData.vp{i}.numLines = VP(i, 5);
end

%% Parse each set of line segments for the corresponding VP
for i = 1: vpData.numVP
    numLine = fscanf(fid, '%d ', [1 1]);
    lines = fscanf(fid, '%g %g %g %g %g', [5 numLine]);
    vpData.vp{i}.lines = lines';
end

end