% crude display of vergini sweep

t = textread('l.gev');
n = size(t,2)-2;        
r = t(:,2);           
k=t(:,1);
gev = t(:,3:n+2);
gev = gev + 1*(gev<1e-16); % kick small values into middle
gev = sort(gev, 2); % sort into ascending order     
subplot(2,1,1);
semilogy(k, 1./gev, '-'); % dir down, neu up
axis([min(k) max(k) 1e-8 1e8]);
subplot(2,1,2);
plot(k, r, '+');
