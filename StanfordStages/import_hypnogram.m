function hypnogram = import_hypnogram(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   RLHYPNOGRAM = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   RLHYPNOGRAM = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   RLhypnogram = importfile('41RL_hypnogram.txt', 1, 296);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/11/11 14:29:34

%% Initialize variables.
delimiter = {''};
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format for each line of text:
%   column1: categorical (%C)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%C%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
hypnogram = table;
hypnogram.score = dataArray{:, 1};

