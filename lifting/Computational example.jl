include("make_start_system.jl")


s = rand(-9:10,21)
Id = identity_matrix(QQ,7)

Gamma = hcat( Id + skew_matrix(s), Id-skew_matrix(s))

SNF0 = hcat(identity(1) + skew_matrix(s), identity(1) - skew_matrix(s))
convert(Matrix{Int64}, SNF0)

g = QQ.(rand(-4:4,7,7).//rand(1:10,7,7))
matrix(QQ,Gamma)
g*Gamma