const assert = require('node:assert/strict');
const test = require('node:test');
const astro = require('../ground-astro.js');

const hours = (h, m, s) => (h + m / 60 + s / 3600) * 15;
const dms = (d, m, s) => Math.sign(d || 1) * (Math.abs(d) + m / 60 + s / 3600);
const near = (actual, expected, tolerance, label) =>
  assert.ok(Math.abs(astro.angleDifference(actual, expected)) <= tolerance,
    `${label}: got ${actual}, expected ${expected} ± ${tolerance}`);
const msFromJulianDate = jd => (jd - 2440587.5) * 86400000;

test('Greenwich mean sidereal time matches Meeus examples 12.a and 12.b', () => {
  // 1987 April 10, 0h UT -> 13h 10m 46.3668s
  near(astro.gmstDegrees(Date.UTC(1987, 3, 10, 0, 0, 0)), hours(13, 10, 46.3668), 1e-5, 'GMST 0h');
  // 1987 April 10, 19h 21m 00s UT -> 8h 34m 57.0896s
  near(astro.gmstDegrees(Date.UTC(1987, 3, 10, 19, 21, 0)), hours(8, 34, 57.0896), 1e-5, 'GMST 19:21');
});

test('precession J2000 -> date matches Meeus example 21.b (theta Persei, 2028)', () => {
  const T = (2462088.69 - 2451545.0) / 36525;
  const [ra, dec] = astro.precessJ2000ToDate(hours(2, 44, 12.975), dms(49, 13, 39.90), T);
  near(ra, hours(2, 46, 11.331), 1 / 3600, 'RA');
  near(dec, dms(49, 20, 54.54), 1 / 3600, 'Dec');
});

test('equatorial -> horizon matches Meeus example 13.b (Venus from Washington)', () => {
  // 1987 April 10, 19:21:00 UT, US Naval Observatory.
  // Meeus measures azimuth westward from south: A = 68.0337, h = 15.1249.
  const ms = Date.UTC(1987, 3, 10, 19, 21, 0);
  const latitude = dms(38, 55, 17);
  const longitude = -dms(77, 3, 56);
  const result = astro.equatorialOfDateToAltAz(hours(23, 9, 16.641), -dms(6, 43, 11.61), latitude, longitude, ms);
  near(result.az, 68.0337 + 180, 0.005, 'azimuth');
  near(result.alt, 15.1249, 0.005, 'altitude');
});

test('zenith, north point and east point land where they should', () => {
  const ms = Date.UTC(2026, 10, 9, 23, 22, 0);
  const latitude = 34.0522, longitude = -118.2437;
  const T = astro.centuriesSinceJ2000(ms);
  const lst = astro.localSiderealDegrees(longitude, ms);

  // Zenith: RA = LST and Dec = latitude, in the frame of date.
  const zenith = astro.altAzToRaDec(0, 90, latitude, longitude, ms);
  const [zRa, zDec] = astro.precessJ2000ToDate(zenith.ra, zenith.dec, T);
  near(zRa, lst, 1e-6, 'zenith RA of date');
  near(zDec, latitude, 1e-6, 'zenith Dec of date');

  // North point of the horizon: Dec = 90 - latitude, 12h from the meridian.
  const north = astro.altAzToRaDec(0, 0, latitude, longitude, ms);
  const [nRa, nDec] = astro.precessJ2000ToDate(north.ra, north.dec, T);
  near(nRa, lst + 180, 1e-6, 'north point RA');
  near(nDec, 90 - latitude, 1e-6, 'north point Dec');

  // East point: on the celestial equator, 6h east of the meridian.
  const east = astro.altAzToRaDec(90, 0, latitude, longitude, ms);
  const [eRa, eDec] = astro.precessJ2000ToDate(east.ra, east.dec, T);
  near(eRa, lst + 90, 1e-6, 'east point RA');
  near(eDec, 0, 1e-6, 'east point Dec');
});

test('alt/az <-> RA/Dec round-trips', () => {
  const ms = Date.UTC(2026, 10, 9, 23, 22, 0);
  for(const [az, alt] of [[0, 0], [45, 20], [123.4, 67.8], [271, -9], [359.9, 89.5], [200, 1]]){
    for(const [latitude, longitude] of [[34.05, -118.25], [-33.87, 151.21], [64.1, -21.9], [0, 0]]){
      const eq = astro.altAzToRaDec(az, alt, latitude, longitude, ms);
      const back = astro.raDecToAltAz(eq.ra, eq.dec, latitude, longitude, ms);
      near(back.alt, alt, 1e-9, `alt ${az}/${alt}`);
      // Azimuth is compared on the sky, so it stays meaningful near the zenith.
      near(back.az * Math.cos(alt * astro.RAD), az * Math.cos(alt * astro.RAD), 1e-7, `az ${az}/${alt}`);
    }
  }
});

test('Polaris sits near the north celestial pole for a northern observer', () => {
  const ms = Date.UTC(2026, 10, 9, 23, 22, 0);
  const latitude = 34.0522;
  const polaris = astro.raDecToAltAz(hours(2, 31, 49.09), dms(89, 15, 50.8), latitude, -118.2437, ms);
  assert.ok(Math.abs(polaris.alt - latitude) < 0.8, `Polaris altitude ${polaris.alt}`);
  assert.ok(Math.abs(astro.angleDifference(polaris.az, 0)) < 1.0, `Polaris azimuth ${polaris.az}`);
});

test('the sky turns 360 degrees in one sidereal day', () => {
  const ms = Date.UTC(2026, 10, 9, 23, 22, 0);
  const siderealDay = 86164.0905 * 1000;
  const a = astro.raDecToAltAz(100, 20, 34.05, -118.25, ms);
  const b = astro.raDecToAltAz(100, 20, 34.05, -118.25, ms + siderealDay);
  near(b.az, a.az, 1e-4, 'azimuth after a sidereal day');
  near(b.alt, a.alt, 1e-4, 'altitude after a sidereal day');
});

test('position angle of the up direction is zero on the meridian', () => {
  // Looking due south on the meridian, "up" points at the north celestial pole.
  const ms = Date.UTC(2026, 10, 9, 23, 22, 0);
  const latitude = 34.0522, longitude = -118.2437;
  const m = astro.horizonToDateMatrix(latitude, longitude, ms);
  const center = astro.multiplyMatrixVector(m, astro.altAzToVector(180, 30));
  const up = astro.multiplyMatrixVector(m, astro.tangentUp(180, 30));
  near(astro.positionAngle(center, up), 0, 1e-9, 'meridian');
  // Looking due east, north is to the left of "up": the zenith lies clockwise
  // from celestial north on the sky, i.e. at a negative (westward) angle.
  const centerEast = astro.multiplyMatrixVector(m, astro.altAzToVector(90, 30));
  const upEast = astro.multiplyMatrixVector(m, astro.tangentUp(90, 30));
  const q = astro.positionAngle(centerEast, upEast);
  assert.ok(q < -30 && q > -90, `east position angle ${q}`);
});

test('the up tangent stays defined at the zenith', () => {
  const up = astro.tangentUp(180, 90); // facing south, looking straight up
  near(up[0], 1, 1e-12, 'points north (behind the observer)');
  near(up[2], 0, 1e-12, 'horizontal');
});

test('refraction is about half a degree at the horizon and tiny overhead', () => {
  const horizon = astro.refractionDegrees(0);
  assert.ok(horizon > 0.44 && horizon < 0.6, `horizon refraction ${horizon}`);
  assert.ok(astro.refractionDegrees(45) < 0.02);
  assert.ok(astro.refractionDegrees(90) < 0.001);
  assert.equal(astro.refractionDegrees(-10), 0);
  assert.ok(astro.apparentAltitude(10) > 10);
});

test('Sun position matches Meeus example 25.a (1992 October 13)', () => {
  const sun = astro.sunEquatorialOfDate(msFromJulianDate(2448908.5));
  near(sun.ra, 198.38083, 0.01, 'Sun RA');
  near(sun.dec, -7.78507, 0.01, 'Sun Dec');
});

test('Sun is on the celestial equator at the March 2026 equinox', () => {
  const sun = astro.sunEquatorialOfDate(Date.UTC(2026, 2, 20, 14, 46, 0));
  assert.ok(Math.abs(sun.dec) < 0.02, `Sun declination ${sun.dec}`);
});

test('Sun is up at local noon and down at local midnight in Los Angeles', () => {
  const noon = astro.sunAltAz(34.0522, -118.2437, Date.UTC(2026, 10, 9, 20, 0, 0));
  const midnight = astro.sunAltAz(34.0522, -118.2437, Date.UTC(2026, 10, 10, 8, 0, 0));
  assert.ok(noon.alt > 30 && noon.alt < 45, `noon altitude ${noon.alt}`);
  assert.ok(midnight.alt < -50, `midnight altitude ${midnight.alt}`);
});
