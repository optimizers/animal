problems = {'small', 'small2', 'medium', 'medium2', ...
            'large', 'large2', 'very', 'very2'};

for k = 1 : length(problems)
  problem = problems{k};
  fprintf('processing problem %s\n', problem);
  [A, b] = hb_to_msm(sprintf('%s.hb', problem));

  % scale problem
  n = size(A, 2);
  s = full(sqrt(sum(A.^2, 1)))';
  s(s == 0) = 1;
  S = spdiags(1./s, 0, n, n);
  A = A * S;

  % scaled minimum least-squares solution
  fprintf('  scaled mls...\n');
  x = mls(A, b);
  fid = fopen(sprintf('%s_scaled_mls.txt', problem), 'w');
  fprintf(fid, '%21.15e\n', x);
  fclose(fid);
end
