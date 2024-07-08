# Hyperelastic-arch-with-contact

Experimenting with using penalty formulation to enforce non-penetration contact condition for a simple arch geometry. The material model is compressible Neohookean hyperelasticity and the solution is discretized with a neural network. The variational energy is computed as a function of the neural network parameters and passed into SQP optimization to find a minimum.
