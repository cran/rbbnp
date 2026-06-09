# rbbnp 1.1.0

This release improves how the frequency window for smoothness estimation is chosen.
By default the confidence bands are now narrower, and regression bands are always
finite. The previous behaviour is still available (see the note below).

## Highlights

- **A new default window rule (`methods_get_xi = "snr"`).** The upper cutoff
  frequency is now chosen from the signal-to-noise ratio of the empirical Fourier
  transform. The original rule of Schennach (2020) relies on a worst-case bound that
  selects no frequency at realistic sample sizes, so the package used to fall back to
  a wide window that gave a flat smoothness envelope and very wide bands. The new
  rule always returns a usable window and produces confidence intervals roughly 0.4
  times as wide.

- **More accurate regression bias bounds.** Two corrections are now on by default.
  `noise_floor = "auto"` uses the noise floor appropriate for a general response, and
  `envelope_use_Y = TRUE` fits the smoothness envelope to the joint spectrum of `Y`
  and `X` rather than the marginal spectrum of `X`, which had under-estimated the bias.

- **Finite bands for difficult cases (`integer_r = TRUE`).** When the data do not show
  a clear power-law decay, the fitted slope can fall below the minimum the method
  assumes. It is now raised to that minimum and the amplitude refit, which keeps the
  bias bound finite. This avoids the very wide bands that could otherwise appear for
  very smooth densities or polynomial conditional means. A warning is shown when this
  adjustment is made.

- **Refreshed plots.** Density, regression, and Fourier-transform plots now use a
  cleaner theme with a colorblind-friendly palette and a legend. The Fourier-transform
  plot shows a wider frequency range with the selected window shaded, so you can see
  where the signal gives way to noise; the new `xi_range` and `expand` arguments
  control the displayed range. Custom colors through `fill_ci` and `fill_bias` still work.

## Notes

- To reproduce the previous behaviour, set `methods_get_xi = "Schennach_loose"`,
  `noise_floor = "compact"`, and `envelope_use_Y = FALSE`.

---

# rbbnp 1.0.0

First stable release with a modern S3 interface, performance improvements, and expanded documentation.

## Highlights

- Stable S3 API for density and regression objects (`bbnp_density`, `bbnp_regression`).
- Bias-bound inference with MSE-optimal bandwidths (no undersmoothing), following Schennach (2020).
- Faster bandwidth selection (FFT-based CV) and improved A/r estimation (LP-based envelope).
- Expanded vignettes and pkgdown website.

---

# rbbnp 0.3.0 (development)

This release focuses on modernizing the package interface and improving performance.

## Highlights

- **Modern S3 interface**: `biasBound_density()` and `biasBound_condExpectation()` now return S3 objects
  (`bbnp_density`, `bbnp_regression`) with standard methods: `print()`, `summary()`, `plot()`, `coef()`,
  `confint()`, and `fitted()` (regression).

- **More faithful implementation of Schennach (2020)**:
  - Improved estimation of smoothness parameters (A, r) using a linear-programming formulation.

- **Performance improvements**:
  - Faster bandwidth cross-validation via an FFT-based approach for larger samples.
  - Vectorized conditional variance computation to reduce bottlenecks.

- **Docs & usability**:
  - Expanded vignettes and pkgdown website organization.

## Notes

- Some internal functions were refactored; most user-facing APIs remain backward compatible.
