function q = memd2(x, varargin)
% The simplest usage: q = memed2(x), where x is composed of n channel of 
% images, its shape is [n,h,w]. Moreover, several parameters can be taken 
% as the input, for example, 
tic
global M N N_dim;
% M, N: height and width of input bidimensioanl data
% N_dim: the number of variates, i.e., the number of bidimensional data,
% number of channel
% ndir: the number of projections
% sd,sd2,tol: step criteria
% stop_crtit£º'stop' or 'fix_h', 'stop' means stopping the sifting process
% accroding to step criteria, and 'fix_h' represents the fixed number of
% IMFs
% stp_cnt: it is enable while 'fix_h' is choosen
% seq: low-discrepancy sequences
% nbit: maximal iteration numbers
%

[x, seq, t, ndir, N_dim, M, N, sd, sd2, tol, nbit, MAXITERATIONS, stop_crit, stp_cnt] = set_value(x, nargin, varargin{:});

global dir_vec;
% convert low-discrepancy sequences to projection directions
dir_vec = get_dir(seq, ndir);
tic

r=x; n_imf=1;
q = cell(1,1);
while ~stop_emd(r, ndir) & n_imf <= 4
    % current mode
    m = r;
    [stop_sift, env_mean] = stop_sifting(m,sd,sd2,tol,ndir);
    
    while ~stop_sift
        m = m - env_mean;
        
        [stop_sift,env_mean] = stop_sifting(m,sd,sd2,tol,ndir);
    end
    
    q{1,n_imf} = m;
    
    n_imf = n_imf + 1;
    disp(n_imf)
    toc
    r = r-m;
end
% Stores the residue
q{1,n_imf} = r;

toc
end

%---------------------------------------------------------------------------------------------------
function stp = stop_emd(r, ndir)
global dir_vec;
ner = zeros(ndir,2);
for it=1:ndir
    % Projection of input signal on n-th (out of total ndir) direction vectors
    y = r .* dir_vec{1,it};
    y = sum(y, 3);
    % Calculates the extrema of the projected signal
    % TODO: it can be further optimized 
    max_logic = imregionalmax(y);
    min_logic = imregionalmin(y);
    % total number of extrema
    ner(it, :) = [sum(max_logic(:)), sum(min_logic(:))];
end

% Stops if the all projected signals have less than 3 extrema
stp = any(ner < 3);
end


% get the direction vectors by seq
function [dir_vec] = get_dir(seq, ndir)
% seq: low-discrepancy sequences
% ndir: number of projections
global M N N_dim;
% to save the direction vectors by cell, the size of each cell is M x N x N_dim
dir_vec = cell(1,ndir);
dir_t = zeros(1,N_dim);

for it=1:ndir
    if (N_dim>3) % Multivariate signal (for N_dim ~=3) with hammersley sequence
        % Linear normalisation of hammersley sequence in the range of -1.00 - 1.00
        b=2*seq(1:end,it)-1;
        
        % Find angles corresponding to the normalised sequence
        tht = atan2(sqrt(flipud(cumsum(b(N_dim:-1:2).^2))),b(1:N_dim-1)).';
        % Find coordinates of unit direction vectors on n-sphere
        dir_t(1:N_dim) = [1 cumprod(sin(tht))];
        dir_t(1:N_dim-1) =  cos(tht) .*dir_t(1:N_dim-1);
        
        for i = 1:N_dim
            dir_vec{1,it}(1:M,1:N,i) = dir_t(i);
        end

    elseif N_dim==3
        % Trivariate signal with hammersley sequence
        % Linear normalisation of hammersley sequence in the range of -1.0 - 1.0
        tt = 2*seq(1,it)-1;
        tt((tt>1))=1;
        tt((tt<-1))=-1;
        
        % Normalize angle from 0 - 2*pi
        phirad = seq(2,it)*2*pi;
        st = sqrt(1.0-tt*tt);
        
        dir_vec{1,it}(1:M,1:N,1)=st * cos(phirad);
        dir_vec{1,it}(1:M,1:N,2)=st * sin(phirad);
        dir_vec{1,it}(1:M,1:N,3)=tt;
    
    else  %Bivariate signal
        dir_vec{1,it}(1:M,1:N,1) = cos(2*pi*it/ndir);
        dir_vec{1,it}(1:M,1:N,2) = sin(2*pi*it/ndir);
    end
end

end

%---------------------------------------------------------------------------------------------------
% computes the mean of the envelopes and the mode amplitude estimate
function [env_mean, nem, amp] = envelope_mean(m,ndir) %new
global M N N_dim;
global dir_vec;
% mean, maximal envelope and corresponding minimal envelope
env_mean=zeros(M,N,N_dim);
env_max=zeros(M,N,N_dim);
env_min=zeros(M,N,N_dim);
amp = zeros(M,N);
for it=1:ndir
    % Projection of input signal on nth (out of total ndir) direction vectors
    y = m .* dir_vec{1,it};
    y = sum(y, 3);
    % Calculates the extrema of the projected signal
    max_logic = imregionalmax(y);
    min_logic = imregionalmin(y);
    
    nem = sum(max_logic(:)) + sum(min_logic(:));
    
    max_index = find(max_logic == 1);
    min_index = find(min_logic == 1);
    
    [max_x, max_y] = ind2sub([M N],max_index);
    [min_x, min_y] = ind2sub([M N],min_index);
    
    for i = 1:N_dim
        max_value = m(max_index + (i-1)*M*N);
        min_value = m(min_index + (i-1)*M*N);
        env_max(:,:,i) = gridfit(max_y,max_x,max_value,1:N,1:M);
        env_min(:,:,i) = gridfit(min_y,min_x,min_value,1:N,1:M);
    end
    
    amp = amp + sqrt(sum((env_max - env_min).^2, 3));
    
    env_mean = env_mean + (env_max + env_min)/2;
end
% average to obtain the final envelope
env_mean = env_mean / ndir;
end


%-------------------------------------------------------------------------------
% Stopping criterion
function [stp,env_mean] = stop_sifting(m,sd,sd2,tol,ndir)
global M N N_dim;
try
    [env_mean,nem,amp] = envelope_mean(m,ndir);
    sx = sqrt(sum(env_mean.^2,3));
    if(amp) % something is wrong here
        sx = sx./amp;
    end
    stp = ~((mean(sx > sd) > tol | any(sx > sd2)) & any(nem > 9));
catch
    env_mean = zeros(M,N,N_dim);
    stp = 1;
end
end

% generate Hammersley sequences
function seq = hamm(n,base)
seq = zeros(1,n);
if ( 1 < base )
    seed = 1:1:n;
    base_inv = inv(base);
    while ( any ( seed ~= 0 ) )
        digit = mod (seed(1:n), base);
        seq = seq + digit * base_inv;
        base_inv = base_inv / base;
        seed = floor (seed / base );
    end
else
    temp = 1:1:n;
    seq = (mod(temp,(-base + 1 ))+0.5)/(-base);
end 
end


% handle the input parameters
function [q, seq, t, ndir, N_dim, M, N, sd, sd2, tol, nbit, MAXITERATIONS, stp_crit, stp_cnt] = set_value(q, narg, varargin)

error(nargchk(1,4,narg));
ndir = [];
stp_crit = [];
stp_vec = [];
stp_cnt  = [];
MAXITERATIONS  = [];
sd=[];
sd2=[];
tol=[];
prm= [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149];

% Changes the input vector to double vector
q = double(q);

% Specifies maximum number of channels that can be processed by the code
% Its maximum possible value is 32.
Max_channels = 16;

if(narg==2)
    ndir=varargin{1};
end

if(narg==3)
    if(~isempty(varargin{1}))
        ndir=varargin{1};
    else
        ndir=64;
    end
    stp_crit=varargin{2};
end

if(narg==4 && strcmp(varargin{2},'fix_h'))
    if(isempty(varargin{1}))
        ndir=64;
        stp_crit=varargin{2};
        stp_cnt  = varargin{3};
    else
        ndir=varargin{1};
        stp_crit=varargin{2};
        stp_cnt  = varargin{3};
    end
elseif (narg==4 && strcmp(varargin{2},'stop'))
    if(isempty(varargin{1}))
        ndir=64;
        stp_crit=varargin{2};
        stp_vec=varargin{3};
    else
        ndir=varargin{1};
        stp_crit=varargin{2};
        stp_vec=varargin{3};
    end
elseif (narg==4 && ~xor(strcmp(varargin{2},'fix_h'),strcmp(varargin{2},'stop')))
    Nmsgid = generatemsgid('invalid stop_criteria');
    error(Nmsgid,'stop_criteria should be either fix_h or stop');
end

%%%%%%%%%%%%%% Rescale input signal if required
if (any(size(q)) == 0)
    datamsgid = generatemsgid('emptyDataSet');
    error(datamsgid,'Data set cannot be empty.');
end


%%%%%%%%%%%% Dimension of input signal
N_dim = size(q,3);
if(N_dim < 2 || N_dim > Max_channels)
    error('Function only processes the signal having 2 and 16 channels.');
end

%%%%%%%%%%%% Length of input signal
M = size(q,1);
N = size(q,2);

%%%%%%%%%%%%% Check validity of Input parameters
if ~isempty(ndir) && (~isnumeric(ndir) || ~isscalar(ndir) || any(rem(ndir,1)) || (ndir < 6))
    Nmsgid = generatemsgid('invalid num_dir');
    error(Nmsgid,'num_dir should be an integer greater than or equal to 6.');
end

if ~isempty(stp_crit) && (~ischar(stp_crit) || ~xor(strcmp(stp_crit,'fix_h'),strcmp(stp_crit,'stop')))
    Nmsgid = generatemsgid('invalid stop_criteria');
    error(Nmsgid,'stop_criteria should be either fix_h or stop');
end

if ~isempty(stp_vec) && (~isnumeric(stp_vec) || length(stp_vec)~=3 || ~strcmp(stp_crit,'stop'))
    Nmsgid = generatemsgid('invalid stop_vector');
    error(Nmsgid,'stop_vector should be an array with three elements e.g. default is [0.075 0.75 0.075] ');
end

if ~isempty(stp_cnt) && (~isnumeric(stp_cnt) || ~isscalar(stp_cnt) || any(rem(stp_cnt,1)) || (stp_cnt < 0) || ~strcmp(stp_crit,'fix_h'))
    Nmsgid = generatemsgid('invalid stop_count');
    error(Nmsgid,'stop_count should be a nonnegative integer');
end

if (isempty(ndir))
    ndir=8; % default
end

if (isempty(stp_crit))
    stp_crit='stop'; % default
end

if (isempty(stp_vec))
%     stp_vec=[0.075,0.75,0.075]; % default
    stp_vec = [0.01, 0.1, 0.01];
end

if (isempty(stp_cnt))
    stp_cnt=2; % default
end

if(strcmp(stp_crit,'stop'))
    sd = stp_vec(1);
    sd2 = stp_vec(2);
    tol = stp_vec(3);
end

%%%%%%%%%%%%% Initializations for Hammersley function
base(1) = -ndir;

%%%%%%%%%%%%%% Find the pointset for the given input signal
if(N_dim==3)
    base(2) = 2;
    for it=1:N_dim-1
        seq(it,:) = hamm(ndir,base(it));
    end
elseif N_dim>3
    for iter = 2 : N_dim
        base(iter) = prm(iter-1);
    end
    
    for it=1:N_dim
        seq(it,:) = hamm(ndir,base(it));
    end
else
    seq=[];
end

%%%%%%%%%%%% Define t
t=1:N;

% Counter
nbit=0;
MAXITERATIONS=1000; % default

% tic
end