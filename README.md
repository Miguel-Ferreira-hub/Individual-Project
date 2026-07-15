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
