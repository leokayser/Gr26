quotMinPol = method()
quotMinPol Ideal := m -> quotMinPol(m, degree m)
quotMinPol (Ideal, ZZ) := (m,d) -> (
    S := ring(m);
    kk := coefficientRing(S);
    alpha = (ring m)_0;
    for i from 1 to 10 do (
        ev := map(S/m, kk[t], {sub(alpha, S/m)});
        mu := (entries mingens ker(ev))_0_0;
        if (degree mu)_0 == d then
            return mu;
        if i == 1 then print "Given element is not a generator. Choosing different element.";
        alpha = random(1, S) + random(0, S);
        print alpha;
    );
    error "No generator found!"
)

dehomogenize = method()
dehomogenize Ideal := I -> dehomogenize(I, (ring I)_0)
dehomogenize (Ideal, RingElement) := (I,x) -> (
    S = ring(I);
    kk = coefficientRing(S);
    S' = kk[delete(x,gens S)];
    I' = sub(I,x=>1);
    return sub(I',S');
)