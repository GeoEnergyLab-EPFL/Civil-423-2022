%% Geometry
m=2;
m1=2;
hd=40;
hw=35;
b=8;
%% Graphical solution

figure(1)
h_axes=0:hd/100:hd;
a0_1=b/m1+hd-sqrt((b/m1+hd)^2-h_axes.^2);
a0_2=m1/m*(hw-h_axes).*log(hd./(hd-h_axes));
plot(h_axes,[a0_1;a0_2]);
title('Graphical solution - Intersection of curves')
a.TickLabelInterpreter = 'latex';
legend('Curve 1','Curve 2')

%% Fixed point solution
tolerance = 1.e-10;
iter_max=10;
k=0;
h=0.6*hw;
err=1;
while (k<iter_max) && (err>tolerance)
 k=k+1;
 h0=h;
 h=hw-m/(m1*log(hd/(hd-h0)))*(b/m1+hd-sqrt((b/m1+hd)^2-h0^2));
 err=abs((h-h0)./h);
 disp(['step number ', num2str(k), ', error = ',num2str(err)]);
end
a0=m1/m*(hw-h)*log(hd/(hd-h));
disp(['h = ', num2str(h),' and a0 = ', num2str(a0)]);