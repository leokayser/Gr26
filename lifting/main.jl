Γ = get_SA_points();# 14x7 rational or ComplexF64 matrix
                    # At least two reasonable approaches:
                    # 1. Get random O(7) matrix, take [Id_7|O(7)]
                    # 2. Get random lin.space [L] in Gr(7,Skew(6)), slice

target = to_A6(Γ);  # Compute target in A6 to which we do the homotopy from start point
                    # This includes normalizing into [Id_7,O(7)]
                    # target ∈ Skew(6), potentially normalized using {±1}^n ⋊ S_7
                    # ((P^7)^14\Δ) / GL(7) ≅ O(7) / ({±1}^n ⋊ S_7) bir'l ≈ Skew(7) / something

function do_the_homotopy(target)
    # load start system
    # load polynomial system
    # run HC
    # collect result
    # Meaningful comment
    return

Λ = do_the_homotopy(target) # Load start system and polynomial system, run HC
                            # This should be a 7x15 matrix, or a 7xSkew(6) thingy

represent_solution(Λ)   # Print plücker coords or basis or whatever of L

validate(Λ,Γ)           # Test if slicing X∩Λ equals Γ after coordinate transformation