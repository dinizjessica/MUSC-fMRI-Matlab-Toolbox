function [W0,R] = null_model_und_sign_helper(CM,opts)

bin_swaps = 5;      % default values defined
wei_freq = 1;       % by null_model_und_sign.m
                    % (BCT toolbox)

if isfield(opts,'bin_swaps')
    if isnumeric(opts.bin_swaps), bin_swaps = opts.bin_swaps; 
    else
        error('myApp:argChk','bin_swaps must be a number.');
    end
end
if isfield(opts,'wei_freq')
    if isnumeric(opts.wei_freq) && (opts.wei_freq>0 && opts.wei_freq<=1)
        wei_freq = opts.wei_freq; 
    else
        error('myApp:argChk','wei_freq must be in the range of: 0 < wei_freq <= 1.');
    end
end

[W0,R] = null_model_und_sign(CM,bin_swaps,wei_freq);

end