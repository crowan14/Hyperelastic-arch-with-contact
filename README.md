# Hyperelastic-arch-with-contact

Experimenting with using penalty formulation to enforce non-penetration contact condition for a simple arch geometry. The material model is compressible Neohookean hyperelasticity and the solution is discretized with a neural network. The variational energy is computed as a function of the neural network parameters and passed into SQP optimization to find a minimum.

The energy functional whose minimum corresponds to a solution is 

$$ \Pi = \int \Psi d\Omega - \int t_i u_i dS + \lambda \int I(g(s))^2 ds $$ 

where $\Psi$ is the strain energy density, $\lambda$ is a contact penalty parameter, $g(s)$ is a function measuring the distance to contact along an a priori known contact surface (gap function), and $I(x)$ is an indicator function which satisfies $I(x)=0$ for $x>0$ and $I(x)=x$ for $x\leq 0$. This simply says that there is a penalty to the energy when the structure penetrates the contact surface.
