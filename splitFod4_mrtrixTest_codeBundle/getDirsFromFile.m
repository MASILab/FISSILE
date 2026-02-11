function dirs_m = getDirsFromFile(dirsFile_s, tabFlag)
% Usage: dirs_m = getDirsFromFile(dirsFile_s, tabFlag)
%
% getDirsFromFile reads gradient directions from a text file
% (made by Mathematica, for example) that contains x, y, and z
% coordinates for one direction on each row of the file. The
% function returns the directions in an nDirs x 3 array. 
% If tabFlag = 1, then it is assumed that x, y, and z are 
% separated by tabs (ASCII 9) on each line. Otherwise (e.g., 
% tabFlag = 0), spaces separate the coordinates on each line.

% Preliminaries:
asciiTab = 9;

% Read directions file:
dirsFid = fopen(dirsFile_s);
dirs_m = [];
while 1
	line_s = fgetl(dirsFid);
	if ~isstr(line_s)
		break
	end
	% Parse line into x, y, and z components. Fields are delimited
	% by spaces or a tab:
    if (tabFlag == 1)
        tabIndex_v = find(abs(line_s) == asciiTab);
        % Now get number fields:
        nx = str2num(line_s(1:(tabIndex_v(1)-1)));
        ny = str2num(line_s((tabIndex_v(1)+1):(tabIndex_v(2)-1)));
        nz = str2num(line_s((tabIndex_v(2)+1):length(line_s)));
    else
        spaceIndex_v = findstr(line_s, ' ');
        % Find indices of first space in a series:
        firstLogic_v = find([2, diff(spaceIndex_v)] > 1); % Must take 1st.
        firstSpaceIndex_v = spaceIndex_v(firstLogic_v);
        % Find indices of last space in a series:
        lastLogic_v = find([diff(spaceIndex_v), 2] > 1);  % Must take last.
        lastSpaceIndex_v = spaceIndex_v(lastLogic_v);
        % Now get number fields:
        nx = str2num(line_s(1:(firstSpaceIndex_v(1)-1)));
        ny = str2num(line_s((lastSpaceIndex_v(1)+1):(firstSpaceIndex_v(2)-1)));
        nz = str2num(line_s((lastSpaceIndex_v(2)+1):length(line_s)));
    end
	dirs_m = [dirs_m; nx, ny, nz];
end

% Close file:
fclose(dirsFid);
