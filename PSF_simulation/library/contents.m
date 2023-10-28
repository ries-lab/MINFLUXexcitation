%Electromagnetic field in the focus region of a microscope objective.
%Copyright © Marcel Leutenegger, 2003-2007, École Polytechnique Fédérale de Lausanne (EPFL),
%Laboratoire d'Optique Biomédicale (LOB), BM - Station 17, 1015 Lausanne, Switzerland.
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
%
%
%  effInit
%     Initialize the parameters.
%
%  effField
%     Electric field near the focus of an objective with the
%     chirp z-transform implementation of the Debye integral.
%
%  effFieldEx
%     Electric field near the focus of an objective with the
%     chirp z-transform implementation of the extended Debye
%     integral.
%
%  effIntensity
%     Electromagnetic intensity in the focus region.
%
%  effPhase
%     Phase map of the focus field. Unwraps the phase of the
%     focus field by trying to minimize phase jumps greater
%     than pi between adjacent points.
%
%  effPlotIntensity
%     Plot the intensity.
%
%  effPlotQuadrants
%     Plot the intensity in 7 of 8 quadrants.
%
%  effPlotSections
%     Plot cross-sections (principal planes) through the focus.
%
help eff/contents
