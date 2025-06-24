# Synchrotron-cooling test

Here we test synchrotron cooling and radiative balance for a group of particles. The particles should accelerate up to gamma_rad level and then their radiative losses and Lorentz-force gains should balance.

In order to make the test work, you need to remove the particle movement from the ``core/pic/pushers/boris.c++``.

The code differs from a regular pic-loop by:

- electric current is not deposited
- star's atmosphere does not inject particles at the beginning, we inject a fixed number of electrons close to the left-side of the box
