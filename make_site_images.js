// Render the pictures used by the GitHub Pages front page (index.html).
//
//   npm install playwright && npx playwright install chromium   (once)
//   node make_site_images.js
//
// Writes JPEG screenshots of the interactive labs and images/social_card.png
// (1200 x 630) into images/. The labs are opened from disk; nothing is fetched.
'use strict';
const path = require('path');
const fs = require('fs');
const { chromium } = require('playwright');

const ROOT = __dirname;
const OUT = path.join(ROOT, 'images');
const url = (file) => 'file://' + path.join(ROOT, 'web', file);

async function setLab06(page, cfg) {
  for (const [key, value] of Object.entries(cfg)) {
    if (['pattern', 'map', 'constant'].includes(key)) await page.selectOption('#' + key, value);
    else {
      await page.locator('#' + key).fill(String(value));
      await page.locator('#' + key).dispatchEvent('input');
    }
  }
  await page.waitForTimeout(300);
}

(async () => {
  fs.mkdirSync(OUT, { recursive: true });
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1200, height: 800 }, deviceScaleFactor: 1, colorScheme: 'light' });
  const jpeg = { type: 'jpeg', quality: 82 };

  // Labs 01-05: the top of each page.
  const labs = [
    ['01_3sat_jewel.html', 'lab01_3sat_jewel.jpg'],
    ['02_quasicrystal_lab.html', 'lab02_quasicrystal.jpg'],
    ['03_tiling_lab.html', 'lab03_tiling.jpg'],
    ['04_sweep_results.html', 'lab04_sweep_results.jpg'],
    ['05_complex_julia_sat.html', 'lab05_complex_julia_sat.jpg'],
  ];
  for (const [file, out] of labs) {
    await page.goto(url(file));
    await page.waitForTimeout(1500);
    await page.screenshot({ path: path.join(OUT, out), ...jpeg });
  }

  // Lab 06: the two planes for several configurations.
  const jewels = [
    ['lab06_golden_rings.jpg', {}],
    ['lab06_golden_spiral.jpg', { shear: 0.8 }],
    ['lab06_alpha_halfturn.jpg', { shear: 0, constant: 'alpha' }],
    ['lab06_hat_jewel.jpg', { constant: 'phi', pattern: 'hat169' }],
    ['lab06_square_fold.jpg', { pattern: 'rect', map: 'square' }],
    ['lab06_weierstrass.jpg', { map: 'wpTri', pattern: 'tri' }],
  ];
  await page.goto(url('06_conformal_jewels.html'));
  await page.waitForTimeout(500);
  for (const [out, cfg] of jewels) {
    await setLab06(page, cfg);
    await page.locator('.planes').screenshot({ path: path.join(OUT, out), ...jpeg });
  }

  // Social card: the alpha jewel's output plane beside the title.
  await setLab06(page, { map: 'exp', pattern: 'rect', constant: 'alpha', sectors: 12 });
  const jewel = (await page.locator('#outCanvas').screenshot()).toString('base64');
  const card = await browser.newPage({ viewport: { width: 1200, height: 630 }, deviceScaleFactor: 1 });
  await card.setContent(`<!doctype html><html><body style="margin:0;width:1200px;height:630px;display:grid;
    grid-template-columns:560px 1fr;align-items:center;background:#f8f7f3;color:#23282a;font-family:system-ui,sans-serif">
    <img src="data:image/png;base64,${jewel}" style="width:520px;height:520px;margin-left:40px;border-radius:12px">
    <div style="padding:0 56px 0 24px">
      <div style="color:#87661c;font-size:22px;letter-spacing:.14em;text-transform:uppercase">A mathematical playground</div>
      <div style="font-size:76px;font-weight:650;line-height:1.02;letter-spacing:-.03em;margin:18px 0 22px">Black Swan<br>Jeweled Math</div>
      <div style="font-size:27px;line-height:1.35;color:#59615f">SAT landscapes, golden-ratio rotations, aperiodic tilings and two complex planes, with controls and independent audits.</div>
    </div></body></html>`);
  await card.screenshot({ path: path.join(OUT, 'social_card.png') });
  await browser.close();
  console.log('wrote', fs.readdirSync(OUT).join(', '));
})();
