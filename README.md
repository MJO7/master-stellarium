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

## Ground view (FOV calculator)

The FOV calculator has an optional **Ground view** (button at the top right of the sky, or the "Ground view" sidebar section). It shows the sky as an observer standing on Earth sees it, in the style of Stellarium:

- Ground and a hills/treeline silhouette below the horizon, with compass points (N, NE, E…) on it. Stars, constellation lines, names and figures below the horizon are hidden.
- The sky is placed for your **location** (device location or typed latitude/longitude/elevation) and a **date and time** you can set, run in real time, speed up, reverse or pause.
- **Sun, Moon and planets** are drawn where they are in your sky, with name labels. The Moon shows its phase with the bright limb tilted correctly. The sidebar lists where each one is (or that it has set); click a row to turn toward it.
- **Time buttons on the sky:** slower/reverse, a large Pause/Play button, faster, and Now. <kbd>Space</kbd> also pauses and resumes.
- A status bar shows location, field of view, the azimuth/altitude you are looking at, and the local date/time with its UTC offset.
- Optional true horizon line, and an optional daylight sky.
- Dragging and the arrow keys change azimuth and altitude with the horizon kept level. <kbd>Shift</kbd> + <kbd>←</kbd>/<kbd>→</kbd> adds a deliberate roll ("Level horizon" removes it). Zoom keys and pinch work as usual.
- Camera, OAG and mosaic frames stay true to size and turn with the sky. Leaving Ground view restores the previous view exactly.

The standard view is the default and is unchanged when Ground view is off.

Positions come from `ground-astro.js` (mean sidereal time, IAU 1976 precession, alt/az ↔ J2000, refraction, and a compact Sun/Moon/planet ephemeris after Paul Schlyter). Stars are placed to about 0.01°, agreeing with JPL Horizons to under 20 arcseconds for Los Angeles on 2026-11-09. The Sun, Moon and planets agree with JPL Horizons to under 2 arcminutes at the two epochs in the tests. The Sun and Moon are drawn at true size but never smaller than a few pixels, planets are drawn as dots sized by brightness, and the landscape is illustrative (the true horizon is at 0° altitude).

## Keyboard (FOV calculator)

| Keys | Action |
|---|---|
| ↑ ↓ ← → | Move the view (in Ground view: altitude and azimuth) |
| Shift + → / Shift + ← | Rotate clockwise / counter-clockwise |
| fn + ↑ / fn + ↓ (Page Up / Page Down) | Zoom in / out |
| Space (Ground view) | Pause / resume time |

## Model notes

The pupil-integration engine has a regression check against a Touch-N-Stars/N.I.N.A. reference calculation. The Antlia comparison is **not** a measurement or a product certification: Antlia's numerical angle-dependent curves and effective coating index are unavailable, so the site labels its assumptions and sensitivity checks. Model percentages are relative to each filter's own peak unless explicitly noted.

## Development

Open the HTML files through any static web server, for example `python3 -m http.server`.

```
node --test tests/*.cjs                                   # run the calculator and astronomy tests
node scripts/build-site.mjs                               # build the static pages
```
