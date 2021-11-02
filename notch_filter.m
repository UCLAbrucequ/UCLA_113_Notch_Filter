clear;


%Set parameters

% f_sampling is the frequency at which the system is sampling
% the signal

% f_reject is the target frequency that we want to reject
f_sampling = 9600;
f_reject = 1200;
f_high = 1300;
f_low = 1100;


% Omega = 2.*pi.*f./f_sampling;
Omega_i = 2.*pi.*f_reject./f_sampling;
Omega_H = 2.*pi.*f_high./f_sampling;
Omega_L = 2.*pi.*f_low./f_sampling;


% Search boundaries
Search_high = 1.0000;
Search_low = 0.9000;

%%
% Step1: Compute a that can reject the frequency at f_reject and also
% satisfy the restrictions at f_high and f_low.

%The "best" a = 0.9359;
a = BinarySearch(1/sqrt(2), Search_high, Search_low, Omega_H, Omega_L);


%computeH(a, Omega_H);
computeBoundaries(a, Omega_H, Omega_L);


% Step2: Compute the transfer function
C = (1-sqrt(2).*a+a.^2)/ (2-sqrt(2));
%H = (1-sqrt(2).*a+a.^2) .* (z.^2 - sqrt(2).* z + 1)./((2-sqrt(2)).*(z.^2 - sqrt(2).* a .* z + a .^2));
%z = linspace(1, 512, 512);
%H = C .* (z.^2 - sqrt(2).* z + 1)/(z.^2 - sqrt(2).* a .* z + a.^2);
%plot(z, H);

% The DF2 Implementation
% y[n] - sqrt(2).* a .* y[n-1] + a.^2 .* y[n-2] = C .* (x[n] - sqrt(2).* x[n-1] + x[n-2])

% Impulse train
Impulse = [1, zeros(1, 511)];

% record the response
y = zeros(1, 512);

% intermediate w
w = zeros(1, 512);

% The first two loops
w(1) = Impulse(1);
y(1) = w(1) * 0.9429;
w(2) = Impulse(2) + (1.324) * w(1);
y(2) = w(2) * 0.9429 - 1.33 * w(1);
% The rest of the loops
for i= 3:512
    w(i) = w(i-2) .* (-0.876) + w(i-1) .* (1.324) + Impulse(i);
    y(i) = w(i) .* 0.9429 + w(i-1) .* (-1.33) + w(i-2) .* 0.9429;
end

%y = y/max(y);


% Figure1
x_axis = linspace(1, 512, 512);

stem(x_axis, y);
title("Figure1: Stem Plot")

omega_axis = linspace(0, pi);


H = zeros(1, length(omega_axis));
for i = 1:length(omega_axis)
    H(i) = (1-sqrt(2).*a+a.^2) .* (exp(j.*2.*omega_axis(i)) - sqrt(2).* exp(j.*omega_axis(i)) + 1)./((2-sqrt(2)).*(exp(j.*2.*omega_axis(i)) - sqrt(2).* a .* exp(j.*omega_axis(i)) + a .^2));
end
figure();

omega_axis = omega_axis .* f_sampling /(2 .* pi);
plot(omega_axis, abs(H));
title("Figure2: Frequency Response");


sinVec = zeros(1, 10000);
for i = 1:10000
    sinVec(i) = sqrt(2).* sin(2.* pi .* 1200 .* i/9600);
end

y_conv = conv(sinVec, y);

I = 0;
for i = 256:10256
    I = I + (y_conv(i)).^2;
end

I;
SIR = 10000/I;


%%
function a = BinarySearch(threshold, high_bound, low_bound, Omega_H, Omega_L)
    % INPUT: threshold, higher_search_bound, lower_search_bound, Omega_high,
    %       Omega_low
    % OUTPUT: the "best" a

    while (high_bound - low_bound >= 0.0001)
        mid = low_bound + (high_bound - low_bound)./2;
        [H_high, H_low] = computeBoundaries(mid, Omega_H, Omega_L);
        
        if (H_high > threshold & H_low > threshold)
            high_bound = mid - 0.00001;
        elseif (H_high <= threshold & H_low <= threshold)
            low_bound = mid + 0.00001;
        end
    end
    a = (low_bound + high_bound)/2;
    
end

 
function [MagOfH] = computeH(a, Omega)
    % INPUT: a, Omega
    % OUTPUT: The magnitude of H
    
    H = (1-sqrt(2).*a+a.^2) .* (exp(j.*2.*Omega) - sqrt(2).* exp(j.*Omega) + 1)./((2-sqrt(2)).*(exp(j.*2.*Omega) - sqrt(2).* a .* exp(j.*Omega) + a .^2));
    MagOfH = abs(H);
end


function [H_high, H_low] = computeBoundaries(a, Omega_H, Omega_L)
    % INPUT: a, Omega_high, Omega_low
    % OUTPUT: The amplitude of H at f_high and f_low
    
    H_high = computeH(a, Omega_H);
    H_low = computeH(a, Omega_L);
end


