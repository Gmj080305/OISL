% RUN_IN_OCTAVE  Run oisl_kk_txcarrier.m in GNU Octave (tested with Octave 8.4).
%
%   MATLAB allows local functions at the END of a script; Octave needs them to be
%   defined BEFORE they are used. This wrapper writes a temporary copy of the
%   script with the function section moved to the top, then runs it.
%   The analysis itself is unchanged. In MATLAB, run oisl_kk_txcarrier.m directly.
%
%   Usage (from this folder):   octave --no-gui run_in_octave.m
%   Figures: needs a graphics toolkit (qt or gnuplot). For a headless run, set
%   PLOT.save = true in oisl_kk_txcarrier.m and start Octave with xvfb-run.

src_file = 'oisl_kk_txcarrier.m';
marker   = '. FUNCTIONS ===';            % header of the local-function section
src = fileread(src_file);
k = strfind(src, marker);
if ~isempty(k), k = find(src(1:k(1)) == sprintf('\n'), 1, 'last') + 1; end   % start of that line
if isempty(k)
    error('Function-section marker not found in %s; update ''marker'' in run_in_octave.m.', src_file);
end
body = strrep(src(1:k(1)-1), 'clear; clc; close all;', 'close all;');
tmp  = [tempname() '_oisl_kk.m'];
fid  = fopen(tmp, 'w');
fputs(fid, sprintf('1;\n'));       % marks the file as a script, not a function file
fputs(fid, src(k(1):end));
fputs(fid, sprintf('\n'));
fputs(fid, body);
fclose(fid);
fprintf('Running %s via %s\n', src_file, tmp);
source(tmp);
delete(tmp);
