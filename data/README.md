### Description of the data

 * real: runs with the real database
 * signal: fake signals with sampling, hiscore signal injected once. No fudging, capping at 2sigma, maxmassdist = 400 gev
 * fake: fake backgrounds, "vanilla" (e.g. no fudging)
 * fudged: same as fake, except fudged by 0.65
 * frozen: fake signals without sampling, hiscore signal injected once. No fudging, capping at 2sigma, maxmassdist = 400 gev.
 * scaled: fake signals, bgerrs scaled by 0.65, hiscore signal injected once. capping at 2sigma, maxmassdist = 400 gev.
 * double: like signal, except hiscore signal injected twice. 
 * iced: like frozen, except hiscore signal injected twice. 
 * cryo: like frozen, except maxmassdist = 100 gev
