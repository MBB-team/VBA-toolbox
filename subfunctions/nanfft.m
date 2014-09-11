function [Y,N,N2] = nanfft(X,N,DIM);
% NANFFT calculates the Fourier-Transform of X for data with missing values. 
%  NANFFT is the same as FFT but X can contain missing values encoded with NaN.
%  NaN's are skipped, NaN do not result in a NaN output. 
%
%   Y = NANFFT(X)
%   Y = NANFFT(X,N)
%   Y = NANFFT(X,[],DIM)
% 
%   [Y,N] = NANFFT(...)
%       returns the number of valid samples N
%
%
% WARNING: missing values can introduce aliasing - causing unintended results.
%    Moreover, the behavior of bandpass and highpass filters in case of missing values 
%    is not fully understood, and might contain some pitfalls.  
%
% see also: FFT, XCORR, NANCONV, NANFILTER

%	$Id$
%	Copyright (C) 2005,2011 by Alois Schloegl <alois.schloegl@gmail.com>		
%       This function is part of the NaN-toolbox available at 
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/ and 
%	http://octave.svn.sourceforge.net/viewvc/octave/trunk/octave-forge/extra/NaN/inst/

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% warning('NANFFT is experimental. For more details see HELP NANFFT');

NX = isnan(X);
X(NX) = 0; 

if nargin==1,
        Y = fft(X);
        N2 = sum(1-NX); % 
        N = fft(NX);
elseif nargin==2,
        Y = fft(X,N);
        N2 = sum(1-NX); 
        N = fft(NX);
elseif nargin==3,
        Y = fft(X,N,DIM);
        N2 = sum(1-NX,DIM); % 
        N = fft(NX,N,DIM);
end;