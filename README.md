# Purpose
The purpose of this project is to explore and associate lithium-ion battery behaviour to state-of-charge (SoC) through electrochemical impedance spectroscopy (EIS), aiming to establish novel methods of determining SoC. 

The project is split into two parts:

1. Quantification of the diffusion process and its relation to SoC, this is performed via the distribution of relaxation times (DRT), a characteristic distribution of the battery system displaying a peak (response) at present time constants.
2. Proposition of equivalent circuit model (ECM) to work back to SoC, this is informed by the results of the DRT.

# Inverse Problem
To calculate the DRT, impedance data $$Z(\omega)$$, is collected from a battery via low amplitude (0.1C) AC signal through EIS (0.01-10,000Hz). The DRT is found by solving the following inverse problem with the aim of establishing the distribution function $$\gamma(\tau)$$:

$$
Z(\omega)=R_{\mathrm{ohmic}}+R_{\mathrm{pol}}
\int_{0}^{\infty}
\frac{\gamma(\tau)}{1+j\omega\tau}d\tau \to R_{\mathrm{ohmic}}+R_{\mathrm{pol}}
\sum_{k=1}^{n}
\frac{\gamma(\tau_k)}{1+j\omega\tau_k}
$$

Where $$Z(\omega)$$ consists of a real $$Z'(\omega)$$, and imaginary $$Z''(\omega)$$ part:

$$
Z(\omega)=Z'(\omega)+jZ''(\omega)
$$

# DRT and Proposed ECM

<p align="center">
  <img src="Images/ECM.png" alt="Logo" width="500">
</p>

<p align="center">
  <img src="Images/DRT.png" alt="Logo" width="500">
</p>

# Diffusion Process
The diffusion here is quantified by:

1. Diffusion peak magnitude $$\gamma(\tau)$$
2. Diffusion time constant $$\tau$$

<p align="center">
  <img src="Images/Relations.png" alt="Logo" width="500">
</p>


# ECM Fit

$$
\min_{x}
\sum_{f=f_{\mathrm{start}}}^{f_{\mathrm{end}}}
\left(
Z_{\mathrm{ECM}}(x)-Z_{\mathrm{EIS}}
\right)
$$

Where $$Z_ecm$$:

$$
Z_{\mathrm{ECM}}
=
R_{\infty}
+
\sum_{i=1}^{4}
\frac{R_i}{1-j\omega R_i C_i}
-
\frac{
\displaystyle \frac{j}{\omega C_5}\left(R_5+\frac{\sigma(1-j)}{\sqrt{\omega}}\right)
}{
\displaystyle \left(R_5+\frac{\sigma(1-j)}{\sqrt{\omega}}\right)-\frac{j}{\omega C_5}
}
$$

$$
x=
\left[
R_{\infty},
R_{1},
R_{2},
R_{3},
R_{4},
R_{5},
C_{1},
C_{2},
C_{3},
C_{4},
C_{5},
\sigma
\right]
$$
