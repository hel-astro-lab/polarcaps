# Compton-cooling test

Here we test Compton cooling and radiative balance for a group of particles. The particles should accelerate up to gamma_rad level and then their radiative losses and Lorentz-force gains should balance.

In order to make the test work, you need to remove the particle movement from the ``core/pic/pushers/boris.c++``.

You should also set ``bool no_photon_update   = true;`` in ``core/qed/interactions/compton.h``.

To match the analytical estimates better it is also advised to calculate the scattering interactions in the Thomson limit by changing the line with ``float s0 = 0.0`` to ``float s0 = 1.0`` and commenting out the if-else block following that line in ``core/qed/interactions/compton.c++``.

The code differs from a regular pic-loop by:

- electric current is not deposited
- star's atmosphere does not inject particles at the beginning, we inject a fixed number of electrons close to the left-side of the box
