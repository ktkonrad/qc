function [p] = empirical_ambiguity_prob(alpha, N, n, ratio, debug)
    % generate an rpw and count the number of times the sign changes
    % between points spaced alpha apart
    % computes probability 2*p(f(x+a/2,y) < 0, f(x,y) > 0, f(x+a,y) > 0)
    % use rpw sampled with ppw_high points per wavelength
    % compute p from N randomly chosen points
    ppw_low = 2*pi/alpha;
    if nargin < 5
        debug = 0;
    end
    if nargin < 4
        ratio = 20;
    end
    if nargin < 3
        n = 50;
    end
    if nargin < 2
        N = 1000000;
    end
    
    ppw_high = ratio * ppw_low;
    
    rpw = rpw2dsample(n*ratio, ppw_high);
    sign_change_count = 0;
    
    if debug
        figure; imagesc(rpw>0); hold on;
    end
    
    for i=1:N
        x = randi((n-1)*ratio);
        y = randi(n*ratio);
        if sign(rpw(y,x)) == sign(rpw(y,round(x+ratio))) && ...
            sign(rpw(y,x)) ~= sign(rpw(y,round(x+ratio/2)))
            sign_change_count = sign_change_count + 1;
            if debug
                plot(x,y, '+r');
                plot(x+ratio,y, '*r');
                x
                y
            end
        elseif debug
            plot(x,y, '+');
            plot(x+ratio,y, '*');
        end
    end
    
    p = 2 * sign_change_count / N;
end
        