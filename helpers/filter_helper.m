function x_f = filter_helper(x, opts)

% x_f = filter_helper(x, opts)
%
% Input:    x,             vector of signals to be normalized by the filter
%
%           opts,          Matlab structure
%                               opts.f_l        low pass frequency cutoff
%                                               in Hz - required
%                                                   (0 < opts.f_l < 1)
%
%                               opts.f_h        high pass frequency cutoff
%                                               in Hz - required
%                                                   (0 < opts.f_h < 1)
%
%                               opts.f_s        sampling frequency in Hz -
%                                               optional
%                                                   (opts.f_s > 0)
%
%                               opts.n          number of FIR coefficients
%                                                   (default 10)
%
%                               opts.type       filter type - required
%                                                   'fir' (default)
%                                                   'butter'
%
% Output:   x_f,            filtered vector with dimension n (where n is
%                           equal to the number of nodes of the vector of
%                           signals. In other words, the number of rows of
%                           the vector of signals).
%
%
% Example:  x = rand(1,225);
%           opts.f_l = 0.01;
%           opts.n = 48;
%           opts.f_h = 0.1;
%           opts.type = 'butter';
%           y = filter_helper(x, opts);
% 

PLOT=false;

if ~exist('opts','var'), error('myApp:argChk', 'Parameter "opts" (MatLab struct) is required.');
elseif ~isstruct(opts), error('myApp:argChk', 'Parameter "opts" must be a MatLab struct.');
end

if  (~isfield(opts,'f_l')) || (~isfield(opts,'f_h')) || (~isfield(opts,'f_l') && ~isfield(opts,'f_h')),
    error('myApp:argChk', 'Low and High Frequency are required parameters.');
elseif ~isa(opts.f_l,'numeric') || ~isa(opts.f_h,'numeric') || (~isa(opts.f_l,'numeric') && ~isa(opts.f_h, 'numeric'))
    error('myApp:argChk', 'Low and High Frequency must be numbers.');
elseif opts.f_l<0 || opts.f_h<0
    error('myApp:argChk', 'Low and High Frequency must be greater than zero.');
elseif (opts.f_l<0 || opts.f_l>1) || (opts.f_h<0 || opts.f_h>1)
    error('myApp:argChk', 'Low and High Frequency must be a number between zero and one.');
elseif (opts.f_l > opts.f_h) && (opts.f_h < opts.f_l)
    error('myApp:argChk', 'Low and High Frequency are switched.');
elseif(opts.f_l == opts.f_h)
    error('myApp:argChk', 'Low and High Frequency are the same.');
end

if ~isfield(opts,'n') || ~isa(opts.n, 'numeric') || opts.n<0
    opts.n = 10;
end

if isfield(opts,'f_s') 
    if ~isa(opts.f_s, 'numeric') || opts.f_s<0
        error('myApp:argChk', 'Sample Frequency must be a valid number.');
    end
    f_nyquist = opts.f_s/2;
else 
    f_nyquist = 1;
end

Wn=[opts.f_l/f_nyquist, opts.f_h/f_nyquist];

if  strcmpi('fir',opts.type),
     b = fir1(opts.n,Wn);
     a = 1;
elseif strcmpi(opts.type,'butter')
     [b,a] = butter(opts.n, Wn);
else
    error('myApp:argChk', 'Invalid type.');
end;


if ~exist('x','var'), error('myApp:argChk', 'Parameter "x" (vector) is required.');
elseif ~isnumeric(x) || isempty(x), error('myApp:argChk', 'Parameter "x" must be a valid vector.');
end

m=mean(x);
xm=x-m;

% filt(b,a,x)
% b: numerator
% a: denominator
% x: input data
x_f = filtfilt(b,a,xm) + m;

if PLOT,
    
   figure;
   hold on;
   plot(1:length(x),x,'k');
   plot(1:length(x_f),x_f,'r','lineWidth',2);
    
end

end
