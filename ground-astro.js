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

  return {
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
