# Finite-Element-Analysis
Finite element analysis using deal.II library of c++</br>

1. Consider the following differential equation of elastostatics, in strong form:</br>
Find u satisfying(E A u,x),x+f A= 0,in (0, L),for the following sets of boundary conditions and forcing function ( ̄fandˆfare constants):</br>
(i)u(0) =g1,u(L) =g2,f= ̄f x,</br>
(ii)u(0) =g1,EAu,x=hatx=L,f= ̄f x,</br>
(iii)u(0) =g1,u(L) =g2,f=ˆf x(L−x)</br>
(iv)u(0) =g1,EAu,x=hatx=L,f=ˆf x(L−x)whereE= 1011Pa,A= 10−4m2, ̄f= 1011Nm−4,ˆf= 1012Nm−5,L= 0.1 m,g1= 0,g2= 0.001 m, andh= 106N.</br>

2. Consider (again) the following differential equation of elastostatics, in strong form:</br>
Find u satisfying(E A u,x),x+f A= 0,in (0, L),for the following sets of boundary conditions and forcing function</br> 
( ̄fandˆfare constants):</br>
(i)u(0) =g1,u(L) =g2,f= ̄f x,</br>
(ii)u(0) =g1,EAu,x=hatx=L,f= ̄f x,</br>
(iii)u(0) =g1,u(L) =g2,f=ˆf x(L−x)</br>
(iv)u(0) =g1,EAu,x=hatx=L,f=ˆf x(L−x)</br>
whereE= 1011Pa,A= 10−4m2, ̄f= 1011Nm−4,ˆf= 1012Nm−5,L= 0.1 m,g1= 0,g2= 0.001 m, andh= 106N.</br>
</br>
3. Solve the steady state problem of heat conduction</br>
PDE−∇·j=f</br>
Constitutive relationj=−κ∇u</br>
Neumann b.c.−j·ˆn=hon∂Ωj</br>
Dirichlet b.c.u=gon∂Ωu</br>
with the following boundary conditions using the specified meshes and linear basis functions.</br>  
Use  ̄κ= 385watt.m−1K−1, whereκij=  ̄κδij.  </br>
Assumej= 0 watt.m−2on all edges/surfaces where no temperature/fluxconditions are specified.</br>

  3.1.  (2D Quadrilateral Mesh):  (x∈[0,0.03],y∈[0,0.08], use a 15 x 40 element mesh.)</br>
  Forcing function:f(x,y) = 0 watt.m−3u(x) = 300(1 +c0x) K alongy= 0 m(bottom nodeset)andu(x) = 310(1 + ˆc0x2) K alongy= 0.08 m(top nodeset)</br>
  wherec0=13K.m−1, ˆc0= 8K.m−2.</br>
</br>
  3.2.  (2D Quadrilateral Mesh):  (x∈[0,0.03],y∈[0,0.08], use a 15 x 40 element mesh.)</br>
  Note:You will calculate theL2norm of the error for this problem (the exact solution is given).</br>
  Forcing function:f(x,y) =−10000 watt.m−3u(x) = 100 + (f/4 ̄κ)x2K alongy= 0 m(bottom nodeset);u(x) = 100 + (f/(4 ̄κ))(x2+ 0.0064) K alongy= 0.08 m(top nodeset);u(y) = 100 + (f/(4 ̄κ))y2K alongx= 0 m(left nodeset);</br>
  u(y) = 100 + (f/(4 ̄κ))(y2+ 0.0009) K alongx= 0.03 m(right nodeset).</br>
  For this case the exact solution isu(x,y) = 100 + (f/(4 ̄κ))(x2+y2).  </br>
  (You can observe the convergenceof your finite element solution with mesh refinement - theL2norm for 5 x 13 elements is 4.1078e-06.)</br>
</br>
4. Solve the steady state problem of heat conduction</br>
PDE−∇·j=f</br>
Constitutive relationj=−κ∇u</br>
Neumann b.c.−j·ˆn=hon∂ΩjDirichlet b.c.u=gon∂Ωu</br>
with the following boundary conditions using the specified meshes and linear basis functions.</br>  
Use  ̄κ= 385watt.m−1K−1, whereκij=  ̄κδij.  </br>
Assumej= 0 watt.m−2on all edges/surfaces where no temperature/fluxconditions are specified.</br>
1. (3D Hexahedral Mesh):  (x∈[0,0.04],y∈[0,0.08],z∈[0,0.02], use a 8 x 16 x 4 element mesh)</br>
Forcing function:f(x,y) = 0 watt.m−3u(y,z) = 300(1 +c0(y+z)) K alongx= 0 m(left nodeset)andu(y,z) = 310(1 +c0(y+z)) K alongx= 0.04 m(right nodeset)</br>
wherec0=13K.m−1.</br>
</br>
5. Consider the 3D elastostatics problem.</br>  
Find u such that</br>
PDEσij,j+fi= 0 in Ω</br>
Constitutive relation σij=Cijklkl</br>
Kinematic relationkl=12(∂uk∂xl+∂ul∂xk)</br>
Neumann b.c.σijnj=hion∂ΩhiDirichlet b.c.ui=ugion∂Ωui</br>
Consider a three-dimensional domain defined byx1= [0,1] m;x2= [0,1] m;x3= [0,1] m (i.e.  the unitcube).</br>   
UseE=  2.0e11  Pa  andν=  0.3.   </br>
Assume  tractionhi=  0  N.m−2on  all  surfaces  where  no  otherconditions are specified.</br>  
Use linear basis functions and a 10 x 10 x 10 element mesh for submission.</br>
Apply the following boundary conditions1.h1=h2= 0,h3= 1.0e9∗x1Pa on the facex3= 1 m;</br>
u1=u2=u3= 0 m on the facex3= 0 m.2.u2= 0.5 + (x2−0.5) cos(π/30)−(x3−0.5) sin(π/30)−x2andu3= 0.5 + (x2−0.5) sin(π/30) + (x3−0.5) cos(π/30)−x3on the facex1= 1 m;u1=u2=u3= 0 m on the facex1= 0 m.</br>

6. Consider a three-dimensional domain defined byx1= [0,1] m;x2= [0,1] m;x3= [0,0.1] m.</br>  
Solve the steadystate and transient heat conduction problems with the following boundary conditions and initial conditions.</br>
Use ρ= 3.8151 x 106N.m−2K−1(specific heat per unit volume),κ= 385 watt.m−1K−1, whereκij=κδij.Assumej= 0 watt.m−2on all edges/surfaces where no temperature/flux conditions are specified.  </br>
Use amesh of 20 x 20 x 1 elements.</br>
1.  (Steady  State  problem):   Boundary  conditionsu=  300  K  alongx1=  0  m  andu=  310  K  alongx1= 1m.2. </br> 
(Transient problem):  Boundary conditionsu= 300 K alongx1= 0 m,u= 310 K alongx1= 1 m.Initial conditionsu= 300 K forx1<0.5 m andu= 300 + 20∗(x1−0.5) K forx1≥0.5 m.</br>

7. Three-dimensional domain defined by [0,1]×[0,1]×[0,1] m (i.e.  the unit cube).  </br>
UseE= 2.0e11Pa,ν=  0.3,ρ=  7.6e3  kg/m3. A  dynamic,  linear  elasticity  code  to  solve  the  following  boundaryvalue problems using linear shape functions.  </br>
Assuming traction hi= 0 N.m−2on all surfaces where no otherconditions are specified.</br>  
Use a 8 x 8 x 8 element mesh.Boundary conditions:u1=u2=u3= 0 m on the facex1= 0 m.</br>
Initial conditions:v1=v2=v3= 0 m/s,u2=u3= 0 m,u1= 0.01x1m.</br>
Use parameters:  ∆t= 1e−6,β= 0.25,γ= 0.5 </br>
