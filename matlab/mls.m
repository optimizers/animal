function x = mls3(A, b)
  F = factorize(A, 'cod');
  x = F \ b;
end
