This repository contains various codes for projects I have been working on with my professor.

# Design of IIR full-band and wide-band integrators with equiripple relative magnitude error 

We have designed integrators up to the fourth order. All poles and zeros of the transfer function are real and simple, situated on the negative part of the real axis.
These integrators, regardless of their order, exhibit an approximately linear phase with a slope of -0.5. In essence, they offer short warm-up periods, and the selection among them depends on the maximum acceptable magnitude and phase errors. 

The approach employed in this paper is based on identifying the initial zeros within the amplitude error and utilizes the Taylor series approximation to achieve an equiripple solution.
The resulting filters demonstrate remarkable performance, surpassing other reported designs in terms of both amplitude and phase error. Additionally, as the order of the proposed integrators increases, both amplitude and phase errors decrease.
It's worth highlighting that the proposed third-order integrator outperforms the fourth-order integrators used for comparison.

On average, all proposed filters exhibit a group delay of 0.5 samples, making them favorable for real-time system applications. Among these designs, the proposed fourth-order integrator demonstrates significantly lower magnitude and phase errors compared to known solutions.
For an integrator of the order N there are 2^(2N) solutions with identical equiripple magnitude errors. The selected proposed integrator is the one with the smallest phase error.

# Comparison of the influence of magnitude and phase error on the error of a differentiator response 

We have investigated how relative magnitude error (up to 10%) and phase error (up to 20 degrees) affects the output of the differentiator. Average odstupanje from the ideal output was calculated.

The input signal was presented in Fourier form: x(t) = sum(C(k)*e^(j*k*w0*t)

Similarly, the output was presented as: y(t) = sum(C'(k)*e^(j*k*w0*t), C' = C*abs(Hid)*e^(j*phid)
Hid and phid are ideal magnitude and phase responses of the diferentiator.

The error signal which is added to the ideal magnitude or phase mimics the one we have attained by placing initial zeros and is a modification of the cosine function.

We have concluded that it is more important to have a smaller phase error than a small amplitude error.

# Paralel all-pass diferentiator corrector

Starting from a digital diferentiator with equiripple magnitude error, we have made a phase corrector with paralel all-pass branches. These correctors were made with orders N = 1, 2, 3 
