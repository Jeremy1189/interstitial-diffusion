function R2 = pearson_R2(outputs,targets)
[N,Q] = size(targets);
m = zeros(N,1);
b = zeros(N,1);
r = zeros(N,1);
for i=1:N
  t = targets(i,:);
  y = outputs(i,:);
  ignore = find(isnan(t) | isnan(y));
  t(ignore) = [];
  y(ignore) = [];
  Quse = Q - length(ignore);
%   h = [t' ones(size(t'))];
%   yt = y';
%   rankStatus = warning('off','MATLAB:rankDeficientMatrix');
%   rankRestore = onCleanup(@() warning(rankStatus));
%   theta = h\yt;
%   m(i) = theta(1);
%   b(i) = theta(2);
  yn = y - mean(y);
  tn = t - mean(t);
  sty = std(yn);
  stt = std(tn);
  r(i) = yn*tn'/(Quse - 1);
  if (sty~=0)&&(stt~=0)
    r(i) = r(i)/(sty*stt);
  end
end
R2=r(i);
end