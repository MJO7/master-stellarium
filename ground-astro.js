// Astronomy helpers for the Ground view: sidereal time, precession, and
// conversions between horizon coordinates (azimuth/altitude) and J2000
// equatorial coordinates (the frame the sky surveys are drawn in).
//
// Conventions
//   - Angles are degrees unless a name says otherwise.
//   - Azimuth is measured from north through east (N = 0, E = 90).
//   - Longitude is positive east, latitude positive north.
//   - Time is a UTC epoch in milliseconds (Date.now()).
//   - Horizon vectors are [north, east, up]; equatorial vectors are
//     [cos(dec)cos(ra), cos(dec)sin(ra), sin(dec)].
//
// Accuracy: mean sidereal time plus IAU 1976 precession. Nutation, annual
// aberration and polar motion are ignored, so positions are good to roughly
// 0.01 degrees, which is far below one screen pixel at any wide field of view.
//
// The file is a classic script so the page can load it directly and the Node
// tests can require() it.
(function(root, factory){
  const api = factory();
  if(typeof module === "object" && module && module.exports) module.exports = api;
  else root.GroundAstro = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function(){
  "use strict";

  const RAD = Math.PI / 180;
  const DEG = 180 / Math.PI;
  const ARCSEC = RAD / 3600;

  function normalizeDegrees(value){
    const v = value % 360;
    return v < 0 ? v + 360 : v;
  }

  // Signed difference a - b wrapped to (-180, 180].
  function angleDifference(a, b){
    let d = (a - b) % 360;
    if(d > 180) d -= 360;
    if(d <= -180) d += 360;
    return d;
  }

  function julianDate(ms){
    return ms / 86400000 + 2440587.5;
  }

  function centuriesSinceJ2000(ms){
    return (julianDate(ms) - 2451545.0) / 36525;
  }

  // Greenwich mean sidereal time (Meeus, Astronomical Algorithms, eq. 12.4).
  function gmstDegrees(ms){
    const d = julianDate(ms) - 2451545.0;
    const T = d / 36525;
    return normalizeDegrees(280.46061837 + 360.98564736629 * d + 0.000387933 * T * T - T * T * T / 38710000);
  }

  function localSiderealDegrees(longitude, ms){
    return normalizeDegrees(gmstDegrees(ms) + longitude);
  }

  // --- small 3x3 helpers (row-major arrays of 9 numbers) -------------------
  function multiplyMatrixVector(m, v){
    return [
      m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
      m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
      m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
    ];
  }

  function transposeMatrix(m){
    return [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]];
  }

  function multiplyMatrices(a, b){
    const out = new Array(9);
    for(let r = 0; r < 3; r++){
      for(let c = 0; c < 3; c++){
        out[r * 3 + c] = a[r * 3] * b[c] + a[r * 3 + 1] * b[3 + c] + a[r * 3 + 2] * b[6 + c];
      }
    }
    return out;
  }

  function dot(a, b){ return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function cross(a, b){
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  }
  function normalize(v){
    const n = Math.hypot(v[0], v[1], v[2]) || 1;
    return [v[0] / n, v[1] / n, v[2] / n];
  }

  // --- spherical <-> vector --------------------------------------------------
  function sphericalToVector(lonDeg, latDeg){
    const lon = lonDeg * RAD, lat = latDeg * RAD, c = Math.cos(lat);
    return [c * Math.cos(lon), c * Math.sin(lon), Math.sin(lat)];
  }

  function vectorToSpherical(v){
    const lon = Math.atan2(v[1], v[0]) * DEG;
    const lat = Math.asin(Math.max(-1, Math.min(1, v[2] / (Math.hypot(v[0], v[1], v[2]) || 1)))) * DEG;
    return [normalizeDegrees(lon), lat];
  }

  // Horizon direction for an azimuth and altitude: [north, east, up].
  function altAzToVector(az, alt){ return sphericalToVector(az, alt); }

  // Unit tangent pointing toward increasing altitude ("up" on a level screen).
  // Stays well defined at the zenith, where it points away from the azimuth
  // the observer is facing.
  function tangentUp(az, alt){
    const a = az * RAD, h = alt * RAD;
    return [-Math.sin(h) * Math.cos(a), -Math.sin(h) * Math.sin(a), Math.cos(h)];
  }

  // Unit tangent pointing toward increasing azimuth ("right" when looking out).
  function tangentRight(az){
    const a = az * RAD;
    return [-Math.sin(a), Math.cos(a), 0];
  }

  // --- precession (IAU 1976), J2000 -> mean equator and equinox of date -----
  function precessionMatrix(T){
    const zeta = (2306.2181 * T + 0.30188 * T * T + 0.017998 * T * T * T) * ARCSEC;
    const z = (2306.2181 * T + 1.09468 * T * T + 0.018203 * T * T * T) * ARCSEC;
    const theta = (2004.3109 * T - 0.42665 * T * T - 0.041833 * T * T * T) * ARCSEC;
    const cz = Math.cos(zeta), sz = Math.sin(zeta);
    const cZ = Math.cos(z), sZ = Math.sin(z);
    const ct = Math.cos(theta), st = Math.sin(theta);
    return [
      cz * ct * cZ - sz * sZ, -sz * ct * cZ - cz * sZ, -st * cZ,
      cz * ct * sZ + sz * cZ, -sz * ct * sZ + cz * cZ, -st * sZ,
      cz * st, -sz * st, ct
    ];
  }

  function precessJ2000ToDate(ra, dec, T){
    return vectorToSpherical(multiplyMatrixVector(precessionMatrix(T), sphericalToVector(ra, dec)));
  }

  // --- horizon <-> equatorial ------------------------------------------------
  // Matrix taking a horizon vector [N, E, Up] to the equatorial frame of date.
  function horizonToDateMatrix(latitude, longitude, ms){
    const lst = localSiderealDegrees(longitude, ms) * RAD;
    const phi = latitude * RAD;
    const sp = Math.sin(phi), cp = Math.cos(phi), sl = Math.sin(lst), cl = Math.cos(lst);
    return [
      -sp * cl, -sl, cp * cl,
      -sp * sl, cl, cp * sl,
      cp, 0, sp
    ];
  }

  // Matrix taking a horizon vector [N, E, Up] to J2000 equatorial coordinates.
  // It is orthonormal, so its transpose is the inverse (J2000 -> horizon).
  function horizonToJ2000Matrix(latitude, longitude, ms){
    const toDate = horizonToDateMatrix(latitude, longitude, ms);
    const dateToJ2000 = transposeMatrix(precessionMatrix(centuriesSinceJ2000(ms)));
    return multiplyMatrices(dateToJ2000, toDate);
  }

  function altAzToRaDec(az, alt, latitude, longitude, ms){
    const m = horizonToJ2000Matrix(latitude, longitude, ms);
    const [ra, dec] = vectorToSpherical(multiplyMatrixVector(m, altAzToVector(az, alt)));
    return {ra, dec};
  }

  function raDecToAltAz(ra, dec, latitude, longitude, ms){
    const m = transposeMatrix(horizonToJ2000Matrix(latitude, longitude, ms));
    const [az, alt] = vectorToSpherical(multiplyMatrixVector(m, sphericalToVector(ra, dec)));
    return {az, alt};
  }

  // Same conversion for coordinates already referred to the equinox of date
  // (used for the Sun, and for checking against published worked examples).
  function equatorialOfDateToAltAz(ra, dec, latitude, longitude, ms){
    const m = transposeMatrix(horizonToDateMatrix(latitude, longitude, ms));
    const [az, alt] = vectorToSpherical(multiplyMatrixVector(m, sphericalToVector(ra, dec)));
    return {az, alt};
  }

  // Position angle (north through east) of a tangent direction at a point on
  // the sphere. Both vectors must be in the same equatorial frame.
  function positionAngle(center, tangent){
    const pole = [0, 0, 1];
    const east = cross(pole, center);
    const eastLength = Math.hypot(east[0], east[1], east[2]);
    if(eastLength < 1e-12) return NaN; // looking exactly at a celestial pole
    const e = [east[0] / eastLength, east[1] / eastLength, east[2] / eastLength];
    const n = cross(center, e);
    return Math.atan2(dot(tangent, e), dot(tangent, n)) * DEG;
  }

  // --- atmospheric refraction ------------------------------------------------
  // Saemundsson's formula: refraction (degrees) for a true (airless) altitude,
  // standard pressure and temperature. About 0.48 degrees at the horizon.
  function refractionDegrees(trueAltitude){
    if(!(trueAltitude > -2)) return 0;
    const h = Math.min(90, trueAltitude);
    const arcmin = 1.02 / Math.tan((h + 10.3 / (h + 5.11)) * RAD);
    return Math.max(0, arcmin) / 60;
  }

  function apparentAltitude(trueAltitude){
    return trueAltitude + refractionDegrees(trueAltitude);
  }

  // --- the Sun (Meeus, chapter 25, low-accuracy method, about 0.01 degrees) --
  function sunEquatorialOfDate(ms){
    const T = centuriesSinceJ2000(ms);
    const L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T * T;
    const M = (357.52911 + 35999.05029 * T - 0.0001537 * T * T) * RAD;
    const C = (1.914602 - 0.004817 * T - 0.000014 * T * T) * Math.sin(M)
      + (0.019993 - 0.000101 * T) * Math.sin(2 * M)
      + 0.000289 * Math.sin(3 * M);
    const omega = (125.04 - 1934.136 * T) * RAD;
    const lambda = (L0 + C - 0.00569 - 0.00478 * Math.sin(omega)) * RAD;
    const epsilon = (23.439291 - 0.0130042 * T + 0.00256 * Math.cos(omega)) * RAD;
    const ra = Math.atan2(Math.cos(epsilon) * Math.sin(lambda), Math.cos(lambda)) * DEG;
    const dec = Math.asin(Math.sin(epsilon) * Math.sin(lambda)) * DEG;
    return {ra: normalizeDegrees(ra), dec};
  }

  function sunAltAz(latitude, longitude, ms){
    const sun = sunEquatorialOfDate(ms);
    return equatorialOfDateToAltAz(sun.ra, sun.dec, latitude, longitude, ms);
  }

  // --- Sun, Moon and planets --------------------------------------------------
  // Paul Schlyter's compact method: orbital elements referred to the equinox of
  // date plus the main perturbations of the Moon, Jupiter, Saturn and Uranus.
  // Good to a few arcminutes for 1900-2100, which is plenty for showing where a
  // planet sits above the horizon. Light-time and aberration are ignored.
  const BODY_NAMES = {
    sun: "Sun", moon: "Moon", mercury: "Mercury", venus: "Venus", mars: "Mars",
    jupiter: "Jupiter", saturn: "Saturn", uranus: "Uranus", neptune: "Neptune"
  };
  const PLANET_IDS = ["mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune"];

  function orbitalElements(id, d){
    switch(id){
      case "sun": return {N: 0, i: 0, w: 282.9404 + 4.70935e-5 * d, a: 1, e: 0.016709 - 1.151e-9 * d, M: 356.0470 + 0.9856002585 * d};
      case "moon": return {N: 125.1228 - 0.0529538083 * d, i: 5.1454, w: 318.0634 + 0.1643573223 * d, a: 60.2666, e: 0.054900, M: 115.3654 + 13.0649929509 * d};
      case "mercury": return {N: 48.3313 + 3.24587e-5 * d, i: 7.0047 + 5.00e-8 * d, w: 29.1241 + 1.01444e-5 * d, a: 0.387098, e: 0.205635 + 5.59e-10 * d, M: 168.6562 + 4.0923344368 * d};
      case "venus": return {N: 76.6799 + 2.46590e-5 * d, i: 3.3946 + 2.75e-8 * d, w: 54.8910 + 1.38374e-5 * d, a: 0.723330, e: 0.006773 - 1.302e-9 * d, M: 48.0052 + 1.6021302244 * d};
      case "mars": return {N: 49.5574 + 2.11081e-5 * d, i: 1.8497 - 1.78e-8 * d, w: 286.5016 + 2.92961e-5 * d, a: 1.523688, e: 0.093405 + 2.516e-9 * d, M: 18.6021 + 0.5240207766 * d};
      case "jupiter": return {N: 100.4542 + 2.76854e-5 * d, i: 1.3030 - 1.557e-7 * d, w: 273.8777 + 1.64505e-5 * d, a: 5.20256, e: 0.048498 + 4.469e-9 * d, M: 19.8950 + 0.0830853001 * d};
      case "saturn": return {N: 113.6634 + 2.38980e-5 * d, i: 2.4886 - 1.081e-7 * d, w: 339.3939 + 2.97661e-5 * d, a: 9.55475, e: 0.055546 - 9.499e-9 * d, M: 316.9670 + 0.0334442282 * d};
      case "uranus": return {N: 74.0005 + 1.3978e-5 * d, i: 0.7733 + 1.9e-8 * d, w: 96.6612 + 3.0565e-5 * d, a: 19.18171 - 1.55e-8 * d, e: 0.047318 + 7.45e-9 * d, M: 142.5905 + 0.011725806 * d};
      case "neptune": return {N: 131.7806 + 3.0173e-5 * d, i: 1.7700 - 2.55e-7 * d, w: 272.8461 - 6.027e-6 * d, a: 30.05826 + 3.313e-8 * d, e: 0.008606 + 2.15e-9 * d, M: 260.2471 + 0.005995147 * d};
      default: return null;
    }
  }

  // Position in the orbit's reference plane (the ecliptic of date).
  function orbitPosition(el){
    const M = normalizeDegrees(el.M) * RAD, e = el.e;
    let E = M + e * Math.sin(M) * (1 + e * Math.cos(M));
    for(let k = 0; k < 8; k++){
      const dE = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
      E -= dE;
      if(Math.abs(dE) < 1e-12) break;
    }
    const xv = el.a * (Math.cos(E) - e);
    const yv = el.a * Math.sqrt(1 - e * e) * Math.sin(E);
    const v = Math.atan2(yv, xv), r = Math.hypot(xv, yv);
    const N = el.N * RAD, i = el.i * RAD, vw = v + el.w * RAD;
    return {
      x: r * (Math.cos(N) * Math.cos(vw) - Math.sin(N) * Math.sin(vw) * Math.cos(i)),
      y: r * (Math.sin(N) * Math.cos(vw) + Math.cos(N) * Math.sin(vw) * Math.cos(i)),
      z: r * Math.sin(vw) * Math.sin(i),
      r, trueLongitude: v + el.w * RAD
    };
  }

  function addToLonLat(p, dLonDeg, dLatDeg, dR){
    const lon = Math.atan2(p.y, p.x) + dLonDeg * RAD;
    const lat = Math.atan2(p.z, Math.hypot(p.x, p.y)) + dLatDeg * RAD;
    const r = p.r + (dR || 0);
    return {x: r * Math.cos(lat) * Math.cos(lon), y: r * Math.cos(lat) * Math.sin(lon), z: r * Math.sin(lat), r};
  }

  function eclipticToEquatorial(x, y, z, obliquityDeg){
    const e = obliquityDeg * RAD;
    return [x, y * Math.cos(e) - z * Math.sin(e), y * Math.sin(e) + z * Math.cos(e)];
  }

  // Geocentric positions, equator and equinox of date. Distances are in
  // astronomical units, except the Moon's, which is in Earth radii.
  function solarSystemPositions(ms){
    const d = julianDate(ms) - 2451543.5;
    const obliquity = 23.4393 - 3.563e-7 * d;
    const sin = a => Math.sin(a * RAD), cos = a => Math.cos(a * RAD);
    const out = {};

    const sunEl = orbitalElements("sun", d);
    const sun = orbitPosition(sunEl);
    const sunX = sun.r * Math.cos(sun.trueLongitude), sunY = sun.r * Math.sin(sun.trueLongitude);
    const sunLonDeg = normalizeDegrees(sun.trueLongitude * DEG);
    {
      const eq = eclipticToEquatorial(sunX, sunY, 0, obliquity);
      const [ra, dec] = vectorToSpherical(eq);
      out.sun = {id: "sun", name: BODY_NAMES.sun, ra, dec, distance: sun.r, magnitude: -26.74, illuminated: 1, phaseAngle: 0,
        angularRadius: Math.asin(0.00465047 / sun.r) * DEG};
    }

    // Moon
    {
      const el = orbitalElements("moon", d);
      const Ms = normalizeDegrees(sunEl.M), Mm = normalizeDegrees(el.M);
      const Ls = Ms + sunEl.w, Lm = Mm + el.w + el.N;
      const D = Lm - Ls, F = Lm - el.N;
      const dLon = -1.274 * sin(Mm - 2 * D) + 0.658 * sin(2 * D) - 0.186 * sin(Ms)
        - 0.059 * sin(2 * Mm - 2 * D) - 0.057 * sin(Mm - 2 * D + Ms) + 0.053 * sin(Mm + 2 * D)
        + 0.046 * sin(2 * D - Ms) + 0.041 * sin(Mm - Ms) - 0.035 * sin(D) - 0.031 * sin(Mm + Ms)
        - 0.015 * sin(2 * F - 2 * D) + 0.011 * sin(Mm - 4 * D);
      const dLat = -0.173 * sin(F - 2 * D) - 0.055 * sin(Mm - F - 2 * D) - 0.046 * sin(Mm + F - 2 * D)
        + 0.033 * sin(F + 2 * D) + 0.017 * sin(2 * Mm + F);
      const dR = -0.58 * cos(Mm - 2 * D) - 0.46 * cos(2 * D);
      const p = addToLonLat(orbitPosition(el), dLon, dLat, dR);
      const eq = eclipticToEquatorial(p.x, p.y, p.z, obliquity);
      const [ra, dec] = vectorToSpherical(eq);
      const moonLon = Math.atan2(p.y, p.x) * DEG, moonLat = Math.atan2(p.z, Math.hypot(p.x, p.y)) * DEG;
      const elongation = Math.acos(Math.max(-1, Math.min(1, cos(sunLonDeg - moonLon) * cos(moonLat)))) * DEG;
      const phaseAngle = 180 - elongation;
      out.moon = {id: "moon", name: BODY_NAMES.moon, ra, dec, distance: p.r, phaseAngle, elongation,
        illuminated: (1 + cos(phaseAngle)) / 2,
        magnitude: -12.7 + 0.026 * Math.abs(phaseAngle) + 4e-9 * Math.pow(phaseAngle, 4),
        angularRadius: Math.asin(0.2725 / p.r) * DEG};
    }

    const Mj = normalizeDegrees(orbitalElements("jupiter", d).M);
    const Msat = normalizeDegrees(orbitalElements("saturn", d).M);
    const Mu = normalizeDegrees(orbitalElements("uranus", d).M);
    for(const id of PLANET_IDS){
      let p = orbitPosition(orbitalElements(id, d));
      if(id === "jupiter"){
        p = addToLonLat(p, -0.332 * sin(2 * Mj - 5 * Msat - 67.6) - 0.056 * sin(2 * Mj - 2 * Msat + 21)
          + 0.042 * sin(3 * Mj - 5 * Msat + 21) - 0.036 * sin(Mj - 2 * Msat) + 0.022 * cos(Mj - Msat)
          + 0.023 * sin(2 * Mj - 3 * Msat + 52) - 0.016 * sin(Mj - 5 * Msat - 69), 0);
      }else if(id === "saturn"){
        p = addToLonLat(p, 0.812 * sin(2 * Mj - 5 * Msat - 67.6) - 0.229 * cos(2 * Mj - 4 * Msat - 2)
          + 0.119 * sin(Mj - 2 * Msat - 3) + 0.046 * sin(2 * Mj - 6 * Msat - 69) + 0.014 * sin(Mj - 3 * Msat + 32),
          -0.020 * cos(2 * Mj - 4 * Msat - 2) + 0.018 * sin(2 * Mj - 6 * Msat - 49));
      }else if(id === "uranus"){
        p = addToLonLat(p, 0.040 * sin(Msat - 2 * Mu + 6) + 0.035 * sin(Msat - 3 * Mu + 33) - 0.015 * sin(Mj - Mu + 20), 0);
      }
      const gx = p.x + sunX, gy = p.y + sunY, gz = p.z;
      const R = Math.hypot(gx, gy, gz), r = p.r;
      const eq = eclipticToEquatorial(gx, gy, gz, obliquity);
      const [ra, dec] = vectorToSpherical(eq);
      const phaseAngle = Math.acos(Math.max(-1, Math.min(1, (r * r + R * R - sun.r * sun.r) / (2 * r * R)))) * DEG;
      const base = 5 * Math.log10(r * R);
      let magnitude;
      if(id === "mercury") magnitude = -0.36 + base + 0.027 * phaseAngle + 2.2e-13 * Math.pow(phaseAngle, 6);
      else if(id === "venus") magnitude = -4.34 + base + 0.013 * phaseAngle + 4.2e-7 * Math.pow(phaseAngle, 3);
      else if(id === "mars") magnitude = -1.51 + base + 0.016 * phaseAngle;
      else if(id === "jupiter") magnitude = -9.25 + base + 0.014 * phaseAngle;
      else if(id === "saturn"){
        // Ring brightness depends on how open the rings appear from Earth.
        const lon = Math.atan2(gy, gx) * DEG, lat = Math.atan2(gz, Math.hypot(gx, gy)) * DEG;
        const ir = 28.06, Nr = 169.51 + 3.82e-5 * d;
        const B = Math.asin(sin(lat) * cos(ir) - cos(lat) * sin(ir) * sin(lon - Nr));
        magnitude = -9.0 + base + 0.044 * phaseAngle - 2.6 * Math.sin(Math.abs(B)) + 1.2 * Math.sin(B) * Math.sin(B);
      }
      else if(id === "uranus") magnitude = -7.15 + base + 0.001 * phaseAngle;
      else magnitude = -6.90 + base + 0.001 * phaseAngle;
      out[id] = {id, name: BODY_NAMES[id], ra, dec, distance: R, magnitude, phaseAngle,
        illuminated: (1 + cos(phaseAngle)) / 2, angularRadius: 0};
    }
    return out;
  }

  // Where the Sun, Moon and planets are in an observer's sky. The Moon is
  // corrected for parallax (it shifts by up to a degree across the Earth).
  function solarSystemHorizon(latitude, longitude, ms){
    const bodies = solarSystemPositions(ms);
    const toHorizon = transposeMatrix(horizonToDateMatrix(latitude, longitude, ms));
    const lst = localSiderealDegrees(longitude, ms);
    const list = [];
    for(const id of ["sun", "moon"].concat(PLANET_IDS)){
      const body = bodies[id];
      let direction = sphericalToVector(body.ra, body.dec);
      let angularRadius = body.angularRadius;
      if(id === "moon"){
        const geocentricLatitude = latitude - 0.1924 * Math.sin(2 * latitude * RAD);
        const rho = 0.99833 + 0.00167 * Math.cos(2 * latitude * RAD);
        const site = sphericalToVector(lst, geocentricLatitude);
        const topo = [
          direction[0] * body.distance - site[0] * rho,
          direction[1] * body.distance - site[1] * rho,
          direction[2] * body.distance - site[2] * rho
        ];
        const topoDistance = Math.hypot(topo[0], topo[1], topo[2]);
        direction = [topo[0] / topoDistance, topo[1] / topoDistance, topo[2] / topoDistance];
        angularRadius = Math.asin(0.2725 / topoDistance) * DEG;
      }
      const vector = multiplyMatrixVector(toHorizon, direction);
      const [az, alt] = vectorToSpherical(vector);
      list.push(Object.assign({}, body, {az, alt, vector, angularRadius}));
    }
    return list;
  }

  return {
    solarSystemPositions, solarSystemHorizon,
    RAD, DEG,
    normalizeDegrees, angleDifference,
    julianDate, centuriesSinceJ2000, gmstDegrees, localSiderealDegrees,
    multiplyMatrixVector, transposeMatrix, multiplyMatrices, dot, cross, normalize,
    sphericalToVector, vectorToSpherical, altAzToVector, tangentUp, tangentRight,
    precessionMatrix, precessJ2000ToDate,
    horizonToDateMatrix, horizonToJ2000Matrix,
    altAzToRaDec, raDecToAltAz, equatorialOfDateToAltAz, positionAngle,
    refractionDegrees, apparentAltitude,
    sunEquatorialOfDate, sunAltAz
  };
});
