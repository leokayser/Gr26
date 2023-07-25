loadPackage "Coisotropy"

--S = QQ[x,y,z,w]
--I = ideal(x*z-y^2, y*w-z^2, x*w-y*z)
S = QQ[x_0..x_5]
I = ideal(x_0^2+x_1^2+x_2^2+x_3^2+x_4^2+x_5^2)

coisotropicForm(I,0, Smooth => true)