loadPackage "Coisotropy"

(k,n) = (2,4)
I = Grassmannian(k-1, n-1, CoefficientRing => QQ, Variable => P);

coisotropicForm(I,0, Smooth => true)