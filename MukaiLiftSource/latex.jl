using Pkg; Pkg.activate("MukaiLiftP6")
using MukaiLiftP6
using HomotopyContinuation
using Oscar
using LinearAlgebra

using DelimitedFiles
using JLD

S_start, L_start = make_start();
parametrized_system, l_start, A_rand = make_poly_system(S_start, L_start);

sol = HomotopyContinuation.read_parameters("Solution_notebook.txt")

L_hat = L_start + sum(sol[i]*A_rand[i] for i=1:69)
round.(L_hat, digits = 3)

using Latexify
latexify(round.(L_hat, digits = 3)) |> print
