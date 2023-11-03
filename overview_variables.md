# Variables for continuum limit

## determination potential

- AIC can be used or not, for now only AIC
- etp: the total error is (m84-m16)/2, with m84 and m16 being the corresponding quantiles. We can take the bootstrapsamples of the potential to be the m50 bootstrapsamples, or we can rescale them so they are errortot/sd(m50) as large. This should include the systematical errors in the bootstrapsamples.
- lowlim(pot), omit: The potential and anisotropy are determined in a range from lower to upper limit, with lower >= 1 and upper <= L/2=8. The lower limit is called lowlim, and it is set to 1, which means the first point is not used. The upper limit is set with omit, and omit number of points are not used. for now, lowlim is always 1 and omit is 0 or 1, so the potential/r0 is calculated in the region 2-8 and 2-7.
- lowlim refers to the region for the anisotropy and lowlimpot refers to the region for the potential, they do not have to be the same, but are by default
- c: This is the parameter is used in the definition of r_0, default -1.65
- normal/sideways: We can determine xi_ren and r_0 from the normal or sideways potential, see proceedings for details of determination of xi.

## determination continuum limit at large volume

- fitlim: points within fitlim of r_0(xi=1) are used, for now fitlim=0.3
- xiconst: we can assume xi(beta)= const or xi(beta)=linear
- degree: degree of the polynomial used for extrapolation to the continuum limit
- lowlimfitpot: Taking the ordered list of anisotropies, we start at the nth one. So 1 means including all anisotropies, and 4 means we include the 4-7th anisotropy.
- xisingle: if this is active, we extrapolate in powers of xi, otherwise, we extrapolate in powers of xi^2

## Extrapolation finite volume:

- direct: Determine beta_ren at L=16, then simulate at those beta_ren and extrapolate with P(L=3) and xi(L=16).
- ratio: determine the ratio of P(L=3) and P(L=16), extrpolate this ratio to 0 with xi(L=16). Then, multiply the result with the continuum limit of L=16.
- for beta, we always take the result at L=16.