"Based on a code from Firedrake"
"Plus moving frame -vv*(phi_x-Delta(phi)_x)"
"Half domain"
"BL: continuous Galerkin"
"rectangle domain"

from firedrake import *
from ini_eta import *
from ini_phi import *
from trans import *
from firedrake.petsc import PETSc
import numpy as np
import time



op2.init()
parameters["coffee"]["O2"] = False

# Timing
t00 = time.time()



""" ________________ Parameters ________________ """

epsilon = 0.05           # Small amplitude parameter
ep=epsilon
mu = 0.0025           # Small dispersion parameter

t=-200
dt = 0.005         # Time step
dt1=dt

""" ___________________ Mesh ___________________ """
mesh = Mesh("dom9_h4.msh")   # Load the mesh file

coords = mesh.coordinates


""" ______________ Function Space ______________ """
V = FunctionSpace(mesh, "CG", 1)    # Vector space

""" ___________ Define the functions ___________ """

eta0 = Function(V, name="eta")
phi0 = Function(V, name="phi")
eta1 = Function(V, name="eta_next")
phi1 = Function(V, name="phi_next")
eta2 = Function(V, name="eta_trans")
phi2 = Function(V, name="phi_trans")

q0 = Function(V)
q1 = Function(V)
phi_h = Function(V)
R_half = Function(V, name="eta_xx")
q_h = Function(V, name="-phi_xx")


q = TrialFunction(V)
v = TestFunction(V)

u1 = Function(V)    # eta(n)
u2 = Function(V)    # eta(n)
u3=Function(V)
u4=Function(V)
ut=Function(V)
uxt=Function(V)
d1 = Function(V)    # eta(n)





""" _____________ Initial solution _____________ """



g=.25
qq=(2/9)**(1/6)*g/np.sqrt(ep)
del1=0.001
ww=107/104
M=qq/(1+2*del1+np.sqrt(ww))
m=M*del1
k1 = -(1+.5*np.sqrt(ww))*M-m
k2 = -(.5*np.sqrt(ww))*M-m
k3 = k2+m
k4 = -k3
k5 = -k2
k6 = -k1

# Expression of eta and phi
x = SpatialCoordinate(mesh)
xx= (x[0]-5)*(4.5)**(1/6)*(ep/mu)**(1/2)
yy= x[1]*(4.5)**(1/3)*ep*(1/mu)**(1/2)
t0 = Constant(t)



eta0=initial_eta(xx,yy,u1,u2,u3,u4,ut,uxt,d1,eta0,k1,k2,k3,k4,k5,k6,t,ep,mu)
phi0=initial_phi(xx,yy,u1,d1,phi0,k1,k2,k3,k4,k5,k6,t,ep,mu)






""" _____________ Weak formulations _____________ """#,time=t

Fphi_h = ( v*(phi_h-phi0)/(0.5*dt) + 0.5*mu*inner(grad(v),grad((phi_h-phi0)/(0.5*dt)))
           + v*eta0 + 0.5*epsilon*inner(grad(phi_h),grad(phi_h))*v)*dx
phi_problem_h = NonlinearVariationalProblem(Fphi_h,phi_h)
phi_solver_h = NonlinearVariationalSolver(phi_problem_h)

# followed by a calculation of a half-step solution :math:`q`, performed using a linear solver::

aq_h = v*q*dx
Lq_h = 2.0/3.0*inner(grad(v),grad(phi_h))*dx


q_problem_h = LinearVariationalProblem(aq_h,Lq_h,q_h)
q_solver_h = LinearVariationalSolver(q_problem_h)




#
Feta = ( v*(eta1-eta0)/dt + 0.5*mu*inner(grad(v),grad((eta1-eta0)/dt))
         - 0.5*((1+epsilon*eta0)+(1+epsilon*eta1))*inner(grad(v),grad(phi_h))
         - mu*inner(grad(v),grad(q_h)))*dx
eta_problem = NonlinearVariationalProblem(Feta,eta1)
eta_solver = NonlinearVariationalSolver(eta_problem)

# and finally the second half-step (explicit this time) for the equation of :math:`\phi` is performed and :math:`q` is computed for the updated solution::
#
Fphi = ( v*(phi1-phi_h)/(0.5*dt) + 0.5*mu*inner(grad(v),grad((phi1-phi_h)/(0.5*dt)))
         + v*eta1 + 0.5*epsilon*inner(grad(phi_h),grad(phi_h))*v)*dx

phi_problem = NonlinearVariationalProblem(Fphi,phi1)
phi_solver = NonlinearVariationalSolver(phi_problem)

Lq = 2.0/3.0*inner(grad(v),grad(phi1))*dx
q_problem = LinearVariationalProblem(aq_h,Lq,q1)
q_solver = LinearVariationalSolver(q_problem)

phi1.assign(phi0)
q_solver.solve()

E_data1 = np.zeros(1)
E1 = assemble( (0.5*eta0**2+0.5*(1+epsilon*eta0)*abs(grad(phi0))**2\
                      +mu*(inner(grad(q1),grad(phi0)))+mu*( - 0.75*q1**2))*dx )
E_data1[0]=E1
            
PETSc.Sys.Print(t,E1)  
""" _____________ Time loop _____________ """

output1 = File('data/data_tes/output.pvd')
output1.write(phi0, eta0,phi2, eta2,  time=t)


# We are now ready to enter the main time iteration loop::
t1=t
step=int(0)
T = int(200)

ddx=1/12
trans_gap=2
    
while t < t1+T:        
      t += dt      
      # PETSc.Sys.Print(t)
      phi_solver_h.solve()
      q_solver_h.solve()
      eta_solver.solve()
      phi_solver.solve()
      q_solver.solve()      
      
      eta0.assign(eta1)
      phi0.assign(phi1)
      step +=int(1)
      if step % 50 == 0:      
           output1.write(phi0, eta0,phi2, eta2,  time=t)
      if step % 400== 0:       
        phi2=trans(coords,phi0,phi2,ddx,trans_gap)
        eta2=trans(coords,eta0,eta2,ddx,trans_gap)
        output1.write(phi0, eta0,phi2, eta2,  time=t)
        eta0.assign(eta2)
        phi0.assign(phi2)        
        PETSc.Sys.Print(t,E1)  
      
      E1 = assemble( (0.5*eta0**2+0.5*(1+epsilon*eta0)*abs(grad(phi0))**2\
                      +mu*(inner(grad(q1),grad(phi0)))+mu*( - 0.75*q1**2))*dx )
      E_data1=np.r_[E_data1,[E1]] 

      np.savetxt('data/data_tes/energy1.csv', E_data1)

    

#Edata.close()
print(time.time() - t00)     # Print computational time (s)
