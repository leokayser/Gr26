module MukaiLiftP6

# Authors: Barbara Betti and Leonie Kayser
# Data: Jun 13, 2024
# Short description: This code accompanies the paper
# [BK24] Barbara Betti and Leonie Kayser, "Mukai lifting of self-dual points in P6", arXiv: https://arxiv.org/abs/2406.02734

using Oscar
using LinearAlgebra
using HomotopyContinuation
using DelimitedFiles
using JLD


export leadexp,
plückercoordinates,
oscar_to_HC_Q,
toric_degen_poly_HC,
poly_to_func,
cayley,
vector_to_skew,
skew_to_vector,
skew_to_SNF,
certify_selfdual,
orthogonal_normal_form,
make_start,
make_poly_system

include("Utils.jl")
include("Lifting.jl")

end # module MukaiLiftP6
