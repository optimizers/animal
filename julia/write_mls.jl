using LinearAlgebra, SparseArrays, DelimitedFiles
using Krylov, Quadmath, HarwellRutherfordBoeing

function scale_ls!(A)
  n, m = size(A)
  s = ones(m)
  for j = 1 : m
    i = A.colptr[j]
    k = A.colptr[j+1] - 1
    nj = (i ≤ k) ? norm(A.nzval[i:k]) : 0.0
    if nj > 0.0
      A.nzval[i:k] ./= nj
      s[j] = nj
    end
  end
  return s
end

for pb in ("small", "small2", "medium", "medium2", "large", "large2", "very", "very2")
  print("$pb: ")
  A = RutherfordBoeingData("../rb/$(pb).rb").data
  b = HarwellBoeingMatrix("../hb/$(pb).hb").rhs[:]
  s = scale_ls!(A)
  A = Float128.(A)
  b = Float128.(b)
  M = A' * A
  rhs = A' * b
  x, stats = minres_qlp(M, rhs, rtol=Float128(1e-20), atol=Float128(1e-25))
  println("✓")
  x = Float64.(x)
  open("../mls/krylov/$(pb)_scaled_mls.txt", "w") do io
    writedlm(io, x)
  end
end
