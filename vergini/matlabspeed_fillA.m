function A = matlabspeed_fillA(a, flag)

for i=1:length(flag)
  if flag(i)>0.5
    A(:,i) = sin(a(:,i)) + cos(a(:,i));
  else
    A(:,i) = sin(a(:,i)/2) + cos(a(:,i)/3);
  end
end
