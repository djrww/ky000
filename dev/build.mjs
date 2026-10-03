// dev/build.mjs — 純 Node 建置驅動（沙箱／CI 用；不需要 moon CLI）
//
// 使用官方編譯器 wasm（@moonbit/moonc-worker）與預編 std .mi/.core（tutuca payload），
// 依序編譯 src/ 下各套件（以 moon.pkg.json 之 import 決定順序），最後連結 JS 目標並執行。
//
//   node dev/build.mjs            # 建置 + 執行 main（義務自證總表）
//   node dev/build.mjs --build    # 只建置
//
// 注意：此檔僅為開發工具；MoonBit 程式庫本身（src/）僅依賴語言內建 std，無第三方依賴。
import fs from "node:fs";
import path from "node:path";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const repo = path.resolve(here, "..");
const require = createRequire(import.meta.url);

let moonc;
try {
  moonc = require("@moonbit/moonc-worker");
} catch (e) {
  console.error("找不到 @moonbit/moonc-worker：請先執行 `node dev/setup.mjs`（需要 npm 網路）。");
  console.error("或使用官方工具鏈：`moon build && moon run src/main`。");
  process.exit(2);
}

const root = path.join(repo, "dev", ".payload", "playground");
if (!fs.existsSync(path.join(root, "manifest.json"))) {
  console.error("缺少預編 std payload：請先執行 `node dev/setup.mjs`。");
  process.exit(2);
}
const manifest = JSON.parse(fs.readFileSync(path.join(root, "manifest.json"), "utf8"));
const target = process.env.MBT_TARGET ?? "js";
const m = manifest.targets[target];
const base = path.join(root, "fs", target);
const load = (p) => [p, new Uint8Array(fs.readFileSync(path.join(base, p)))];
const stdMi = m.std.map(load);
const libMi = m.lib.map(load);
const stdCores = m.linkOrder.map((p) => new Uint8Array(fs.readFileSync(path.join(base, p))));

// ---- 讀取模組設定，依 import 排序套件 ----
const mod = JSON.parse(fs.readFileSync(path.join(repo, "moon.mod.json"), "utf8"));
const moduleName = mod.name;
const srcRoot = path.join(repo, mod.source ?? ".");
const pkgs = [];
for (const dir of fs.readdirSync(srcRoot)) {
  const full = path.join(srcRoot, dir);
  if (!fs.statSync(full).isDirectory()) continue;
  const pkgFile = path.join(full, "moon.pkg.json");
  if (!fs.existsSync(pkgFile)) continue;
  const cfg = JSON.parse(fs.readFileSync(pkgFile, "utf8"));
  pkgs.push({
    dir: full,
    name: `${moduleName}/${dir}`,
    deps: cfg.import ?? [],
    isMain: cfg["is-main"] === true || cfg.is_main === true,
  });
}
const ordered = [];
const seen = new Set();
const visit = (p) => {
  if (seen.has(p.name)) return;
  seen.add(p.name);
  for (const d of p.deps) {
    const q = pkgs.find((x) => x.name === d);
    if (q) visit(q);
  }
  ordered.push(p);
};
for (const p of pkgs.filter((p) => p.isMain).concat(pkgs.filter((p) => !p.isMain))) visit(p);

// ---- 逐套件編譯 ----
const t0 = Date.now();
const mi = [];
const cores = [];
for (const p of ordered) {
  const files = fs
    .readdirSync(p.dir)
    .filter((f) => f.endsWith(".mbt"))
    .sort()
    .map((f) => [f, fs.readFileSync(path.join(p.dir, f), "utf8")]);
  const bp = moonc.buildPackage({
    mbtFiles: files,
    miFiles: mi,
    indirectImportMiFiles: libMi,
    stdMiFiles: stdMi,
    target,
    pkg: p.name,
    pkgSources: [`${p.name}:${p.dir}`],
    isMain: p.isMain,
    errorFormat: "human",
    enableValueTracing: false,
    noOpt: false,
  });
  const diags = (bp.diagnostics || []).join("\n").trim();
  if (diags) console.log(diags);
  if (!bp.core) {
    console.log(`\n[建置失敗] ${p.name}（${Date.now() - t0} ms）`);
    process.exit(1);
  }
  mi.push([`${p.name}.mi`, bp.mi]);
  cores.push(bp.core);
}

// ---- 連結 ----
const mains = ordered.filter((p) => p.isMain);
if (mains.length === 0) {
  console.log(`[建置完成] 函式庫套件 ${ordered.length} 個，無 executable 套件（${Date.now() - t0} ms）`);
  process.exit(0);
}
const lk = moonc.linkCore({
  coreFiles: [...cores, ...stdCores],
  main: mains[0].name,
  pkgSources: ordered.map((p) => `${p.name}:${p.dir}`),
  target,
  exportedFunctions: [],
  outputFormat: "wasm",
  testMode: false,
  debug: false,
  noOpt: false,
  sourceMap: false,
  sources: {},
  stopOnMain: false,
});
const outDir = path.join(repo, "dev", ".out");
fs.mkdirSync(outDir, { recursive: true });
const outPath = path.join(outDir, "ky000.mjs");
fs.writeFileSync(outPath, Buffer.from(lk.result));
console.log(`[建置成功 ${Date.now() - t0} ms] ${ordered.map((p) => p.name).join(" -> ")}`);
if (!process.argv.includes("--build")) {
  await import(outPath);
}
