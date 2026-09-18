// Uses the installed Playwright browser, never system Chrome or user profiles.
// Run the browser-automation-guard audit before and after this script.
const {chromium} = require('playwright');
const fs = require('node:fs');
const path = require('node:path');
const {pathToFileURL} = require('node:url');
const assert = require('node:assert/strict');

(async () => {
  const output = path.resolve(process.argv[2] || 'output/visual-qa');
  const qa = path.join(output, 'browser-qa');
  fs.mkdirSync(qa, {recursive:true});
  const server = await chromium.launchServer({headless:true});
  const owner = {pid:server.process().pid, arguments:server.process().spawnargs,
                 task:'neurotransform visual QA', persistent:false};
  fs.writeFileSync(path.join(qa, 'browser-owner.json'), JSON.stringify(owner,null,2));
  let browser;
  const errors=[];
  const checks=[];
  try {
    browser = await chromium.connect(server.wsEndpoint());
    const page = await browser.newPage({viewport:{width:1500,height:1180},deviceScaleFactor:1});
    page.on('pageerror', e=>errors.push(String(e)));
    page.on('console', m=>{if(m.type()==='error')errors.push(m.text())});
    await page.goto(pathToFileURL(path.join(output,'index.html')).href);
    await page.waitForSelector('#canvas-2-0');
    const ids=await page.evaluate(()=>QA.cases.map(c=>c.id));
    for(const id of ids){
      await page.locator(`.family button[data-id="${id}"]`).click();
      if(await page.locator('.unavailable-panel').count()){
        checks.push({id,unavailable:true});
        await page.screenshot({path:path.join(qa,id+'.png'),fullPage:true});
        continue;
      }
      assert.equal(await page.locator('.tile canvas').count(),12);
      const before=await page.locator('#canvas-2-0').evaluate(c=>c.toDataURL());
      await page.locator('#slice-2').focus();
      await page.keyboard.press('ArrowRight');
      const after=await page.locator('#canvas-2-0').evaluate(c=>c.toDataURL());
      assert.notEqual(before,after,`${id}: slice control must change rendered reference`);
      await page.screenshot({path:path.join(qa, id+'.png'),fullPage:true});
      const variants=await page.locator('#scenario option').evaluateAll(xs=>xs.map(x=>x.value));
      for(const variant of variants){
        await page.selectOption('#scenario',variant);
        for(const diagnostic of ['error','coord_error','jacobian'])await page.selectOption('#diagnostic',diagnostic);
      }
      await page.selectOption('#comparison','overlay');
      await page.locator('#opacity').fill('80');
      await page.locator('#opacity').dispatchEvent('input');
      if(id==='ants_affine_warp')await page.screenshot({path:path.join(qa,'wrong-order.png'),fullPage:true});
      const overflow=await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth);
      assert.equal(overflow,false,`${id}: desktop horizontal overflow`);
      checks.push({id,canvases:12,variants:variants.length,slice_changed:true});
    }
    await page.setViewportSize({width:390,height:844});
    await page.selectOption('#mobile-case','fsl_relative');
    await page.screenshot({path:path.join(qa,'mobile.png'),fullPage:true});
    assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false,'mobile overflow');
    assert.deepEqual(errors,[]);
    fs.writeFileSync(path.join(qa,'checks.json'),JSON.stringify({checks,errors,mobile_width:390},null,2));
    console.log(`Verified ${checks.length} cases, all controls, desktop and mobile; no browser errors.`);
    await page.close();
  } finally {
    if(browser)await browser.close();
    await server.close();
    fs.writeFileSync(path.join(qa,'browser-closed.json'),JSON.stringify({pid:owner.pid,closed:true}));
  }
})().catch(e=>{console.error(e);process.exitCode=1});
