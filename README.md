# OISL
# optical inter satellite link
An Optical Inter-Satellite Link (OISL) is essentially the space equivalent of fiber-optic internet—minus the actual cables.
Instead of using traditional radio waves (RF) to beam data back and forth, satellites equipped with OISLs use infrared lasers to talk directly to one another in the vacuum of space.
## A] What typical system model consists of((project model may differ):
#### 1. The Transmitter Model
This models how the light is generated and sent. It includes parameters like:
* **Laser Power ($P_t$):** How bright the laser is.
* **Wavelength ($\lambda$):** Usually in the infrared spectrum (e.g., 1550 nm).
* **Modulation Scheme:** How data is encoded into the light (e.g., turning the laser on and off rapidly, known as On-Off Keying, or using more advanced phase modulation).

#### 2. The Geometric and Orbital Model
Satellites are constantly moving relative to each other at hypersonic speeds. This part of the model calculates:
* The dynamically changing distance ($d$) between satellites.
* The **Doppler shift** (how the movement stretches or compresses the light waves).
* The precise angles the satellites must tilt to maintain a line of sight.

#### 3. The Channel Model (Free Space)
This models what happens to the laser beam as it travels through the vacuum of space. The most critical factor here is **Free Space Path Loss ($L_s$)**, which calculates how much the laser beam naturally spreads out over huge distances:

$$L_s = \left(\frac{\lambda}{4\pi d}\right)^2$$

It also accounts for **pointing jitter**—the tiny mathematical losses caused by the satellite vibrating or slightly missing its target.

#### 4. The Receiver Model
This models how the receiving satellite captures and decodes the light. It includes:
* **Aperture Size:** The diameter of the collecting telescope.
* **Sensor Sensitivity:** How good the photodetector is at spotting faint photons.
* **Noise Sources:** Filtering out background light from the sun, stars, or Earth.
