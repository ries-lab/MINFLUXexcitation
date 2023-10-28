%Introducing Lateral offset of next element, meaning the next modifier in Ei (misalignment).
%
% sys    System data
%  .Da      Aperture diameter [m]
%  .ox      Lateral offset(s) [m]
%  .oy
%
% E      Electric field [V/m]
% r      Radial position [Da/2]
% t      Incidence angle [rad]
% p      Polar angle [rad]

%Copyright ? Marcel Leutenegger, 2003-2007, ?cole Polytechnique F?d?rale de Lausanne (EPFL),
%Laboratoire d'Optique Biom?dicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.
%
%    This library is free software; you can redistribute it and/or modify it under
%    the terms of the GNU Lesser General Public License as published by the Free
%    Software Foundation; version 2.1 of the License.
%
%    This library is distributed in the hope that it will be useful, but WITHOUT ANY
%    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
%    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public License along
%    with this library; if not, write to the Free Software Foundation, Inc.,
%    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
%
function E=offset(sys,E,r,t,p)
% Take the variable, n, into the current workspace.
n=evalin('caller','n') + 1; % Originally it was +1. Execute MATLAB expression in the caller workspace(The workspace of the function, which called the currently executing function).
'offset start1'

if n <= numel(sys.Ei) % numel(A) returns the number of elements, n, in array A
    
    'offset start2'
    % Shifting 'r' value in XY plane. The actual shift is 'ox' in 'r' value
    % and 'oy' in 'p' value.
    r=cis(p,r) - 2/sys.Da*complex(sys.ox(1),sys.oy(1));     % 'cis' function gives Z, Z = r*exp(i*p). 'complex' function z = a + bi. 
    sys.ox=sys.ox(2:end);   
    sys.oy=sys.oy(2:end);   % Taking the pre-defined offset amount, oy, from its second row.
    p=angle(r); % Taking angle component of 'r' as a phase.
    r=abs(r);   % taking absolute component of 'r' as amplitude.
    E=feval(sys.Ei{n},sys,E,r,t,p);
    assignin('caller','sys',sys);
    assignin('caller','n',n);
    'offset end'
    E;
end
end
