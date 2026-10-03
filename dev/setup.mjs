// dev/setup.mjs — 下載固定版本的官方編譯器 wasm 與預編 std 介面，供 dev/build.mjs 使用。
// 只有開發驗證需要；MoonBit 程式庫本身無第三方依賴。
import { execFileSync } from "node:child_process";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const dev = here;
const worker = "0.1.202608195"; // 必須與 payload 的 mooncWorker 相同（同一 moonc 建置）
const payload = "0.23.0";

const pkgJson = { name: "ky000-dev", private: true, type: "module", dependencies: {} };
fs.writeFileSync(path.join(dev, "package.json"), JSON.stringify(pkgJson, null, 2));
console.log("安裝 @moonbit/moonc-worker@" + worker + " 與預編 payload…");
execFileSync("npm", ["install", "--no-fund", "--no-audit", "--silent",
  `@moonbit/moonc-worker@${worker}`, `@marianoguerra/tutuca-playground-payload@${payload}`],
  { cwd: dev, stdio: "inherit" });

const from = path.join(dev, "node_modules", "@marianoguerra", "tutuca-playground-payload", "playground");
const to = path.join(dev, ".payload", "playground");
fs.rmSync(path.join(dev, ".payload"), { recursive: true, force: true });
fs.cpSync(from, to, { recursive: true });
const manifest = JSON.parse(fs.readFileSync(path.join(to, "manifest.json"), "utf8"));
console.log(`完成：payload ${manifest.toolchain}（moonc ${manifest.mooncBuild}）/ worker ${manifest.mooncWorker}`);
console.log("現在可以執行： node dev/build.mjs");
