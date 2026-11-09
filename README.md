# Master Stellarium

Browser-based astronomy calculators for planning imaging sessions. Everything is static HTML, so there is no backend.

**Live site:** https://mjo7.github.io/master-stellarium/ (published with GitHub Pages via `.github/workflows/pages.yml`; Pages must be set to the "GitHub Actions" source in repo settings.)

## Calculators

- **Narrowband filter calculator** (`narrowband-filter-calculator.html`): compares modeled Hα, OIII, and SII transmission through Antlia Standard and Highspeed 3 nm filters. Optical presets include:
  - CarbonStar 150 at f/4, with its 0.95× PRCC at f/3.8, and with the APM 1.5× coma-correcting Barlow at f/6
  - EdgeHD 9.25 at native f/10, or with its dedicated 0.7× reducer at f/7
  - All listed Samyang 135 mm f/2 aperture-ring settings
- **FOV calculator** (`fov-calculator.html`): field-of-view viewer built on a bundled, fixed Aladin Lite release (in `vendor/`) that uses public HiPS survey mirrors.
- `index.html` is a chooser page linking to both calculators.

## Model notes

The pupil-integration engine has a regression check against a Touch-N-Stars/N.I.N.A. reference calculation. The Antlia comparison is **not** a measurement or a product certification: Antlia's numerical angle-dependent curves and effective coating index are unavailable, so the site labels its assumptions and sensitivity checks. Model percentages are relative to each filter's own peak unless explicitly noted.

## Development

Open the HTML files through any static web server, for example `python3 -m http.server`.

```
node --test tests/narrowband-filter-calculator.test.cjs   # run the calculator tests
node scripts/build-site.mjs                               # build the static pages
```
