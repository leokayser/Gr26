include("../Utilities.jl")


R,p,_,_ = plückercoordinates(2,6,QQ);
T,_ = polynomial_ring(QQ, vcat(["t$i" for i=1:8], "u" ) )
#u = gens(T)[9]
p = hom(R,T,gens(T)[1:8]).(p)
w = [-3;-2;-1;0;0;-1;-2;-3;0];

@var x[1:8] u a[1:8,1:15]
pu_HC = [toric_degen_poly_HC(p[i],w,[x;u],gens(T)[9]) for i=1:15]
F = a*pu_HC
C = System(F, variables = x, parameters = [u;a[:]])

println("Defining system done.")

A_rand    = randn(ComplexF64, length(a))

a0_start  = [0;A_rand]
toric_sol = solutions(HomotopyContinuation.solve(C, target_parameters = a0_start))

a1_start  = [1;A_rand]
start_sol = solutions(HomotopyContinuation.solve(C, toric_sol;
    start_parameters = a0_start, target_parameters = a1_start))

# given A
A    = randn(ComplexF64, length(a))
a1_target = [1;A]
X8_cap_A = solutions(HomotopyContinuation.solve(C, start_sol;
    start_parameters = a1_start, target_parameters = a1_target))

HomotopyContinuation.write_parameters("Gr26_start_parameters.txt", a1_start)
HomotopyContinuation.write_solutions("Gr26_start_system.txt", start_sol)


############

   
a1_start = HomotopyContinuation.read_parameters("Gr26_start_parameters.txt")
start_sol = HomotopyContinuation.read_solutions("Gr26_start_system.txt")


########################

# Check that we have obtained points on the grassmannian by evaluating pluecker relations at those points
 

#numerical_plücker = poly_to_fp(ϕ) 
#Z = numerical_plücker.(gr_start_system) # 14 points in P^14

#[A*z for z in Z] #the points in Z lies in the linear section of P^14 as they should


########################

println("Loading start parameters done!")

function uniform_Gr_point(a)
    return randn(length(a[:]))
end

function count_rand_sample()
    targ_par = [1; uniform_Gr_point(a)]
    slice_result = HomotopyContinuation.solve(C, gr_start_system; start_parameters = gr_start_param, target_parameters = targ_par, show_progress=false)
    if nresults(slice_result, only_nonsingular=true) < 14
        return -1
    end
    return nresults(slice_result, only_real=true)
end

N = parse(Int64, ARGS[1])
#N = 10000
counter = Dict(i => 0 for i=-1:14)
for i=1:N
    r = count_rand_sample()
    counter[r] += 1
    if i%1000 == 0
        #file = open("test_count.txt", "w")
        file = open(ARGS[2]*".txt", "w")
        for i in -1:14
            write(file, string(i) * " " * string(counter[i]) * "\n")
        end
        close(file)
        println(string(i) * "/" * string(N) * " done")
    end
end

#file = open("test_count.txt", "w")
file = open(ARGS[2]*".txt", "w")
for i in -1:14
    write(file, string(i) * " " * string(counter[i]) * "\n")
end
close(file)

