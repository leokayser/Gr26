include("../Utilities.jl")


k = 2
m = 6
R,ϕ,vrs,M = plückercoordinates(k,m,QQ);
T, variable = PolynomialRing(QQ, vcat(["t$i" for i=1:8], "u" ) )
u = variable[9]
ι = hom(R,T,gens(T)[1:8])
ϕ = ι.(ϕ)
w = [-3;-2;-1;0;0;-1;-2;-3;0];

#=
terms = [ collect(Oscar.terms(ϕ[j])) for j in eachindex(ϕ) ]
ϕu = [sum([ u^( transpose(w)*(leadexp(terms[j][i],w)-leadexp(ϕ[j],w) )) * terms[j][i]  for i=1:length(terms[j]) ]) for j=1:length(ϕ)]
#l = 8
#Fu = [rand(-100:100,length(ϕu))'*ϕu for i = 1:l]

@var x[1:9] a[1:8,1:15]
vars_HC =  x[1:9] 
ϕu_HC= [oscar_to_HC_Q(ϕu[i], vars_HC) for i=1:length(ϕu)]
=#
ϕu_HC = [toric_degen_poly_HC(ϕ[i],w) for i=1:15]
Fu_HC = a*ϕu_HC #a are the linear forms cutting out a P^6 in P^14

C = System(Fu_HC, variables = vars_HC[1:8], parameters = [vars_HC[9];a[:]])

println("Defining system done.")

###########

#=
targ_par = [0;randn(ComplexF64, length(a[:]))]
torus_result = HomotopyContinuation.solve(C, target_parameters = targ_par )
torus_sol = solutions(torus_result);

targ_par_new = [1;randn(ComplexF64, length(a[:]))]
grass_result = HomotopyContinuation.solve(C, torus_sol; start_parameters = targ_par, target_parameters = targ_par_new)
grass_sol = solutions(grass_result)
HomotopyContinuation.write_parameters("Gr26_start_parameters.txt", targ_par_new)
HomotopyContinuation.write_solutions("Gr26_start_system.txt", grass_sol)
=#

############

   
gr_start_param = HomotopyContinuation.read_parameters("Gr26_start_parameters.txt")
gr_start_system = HomotopyContinuation.read_solutions("Gr26_start_system.txt")

A = reshape(gr_start_param[2:length(gr_start_param)], 8, 15) # This is the linear system whose solution are the Λ in P^14 
B = LinearAlgebra.nullspace(A) # each column is an element of the basis of the space  Λ in P^14. The basis determines an isomorphism with P^6




########################

# Check that we have obtained points on the grassmannian by evaluating pluecker relations at those points
 

numerical_plücker = poly_to_fp(ϕ) 
Z = numerical_plücker.(gr_start_system) # 14 points in P^14

[A*z for z in Z] #the points in Z lies in the linear section of P^14 as they should




########################

println("Loading start parameters done!")

function uniform_Gr_point(a)
    return randn(length(a[:]))
end

#N = parse(Int64, ARGS[1])
N = 100
counter = Dict(2*i => 0 for i=0:7)

for i=1:N
    targ_par = [1; uniform_Gr_point(a)]
    grass_result = HomotopyContinuation.solve(C, gr_start_system; start_parameters = gr_start_param, target_parameters = targ_par)
    real_sols = length(HomotopyContinuation.real_solutions(grass_result))
    counter[real_sols] += 1
end

file = open("test_count.txt", "w")
#file = open(ARGS[2], "w")
for i in 0:7
    write(file, string(counter[2*i]),"\n")
end
close(file)

