function h = VBA_figure (varargin)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% h = VBA_figure (varargin)
% surcharge of the figure function with some styling
%
% IN: optional arguments for the call to figure()
% OUT: handle to the create figure
%
% /////////////////////////////////////////////////////////////////////////

options = [{ ...
    'Color', 'white', ...
    'MenuBar', 'none' ...
    }, varargin];

h = figure (options{:});