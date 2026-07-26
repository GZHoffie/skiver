import { copyFileSync, mkdirSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { spawnSync } from "node:child_process";

const uiDir = dirname(dirname(fileURLToPath(import.meta.url)));
const repoDir = dirname(uiDir);

function run(command, args, cwd) {
  const result = spawnSync(command, args, { cwd, stdio: "inherit", encoding: "utf8" });
  if (result.error) throw result.error;
  if (result.status !== 0) process.exit(result.status ?? 1);
  return result.stdout;
}

const rustc = spawnSync("rustc", ["-vV"], {
  cwd: repoDir,
  encoding: "utf8",
});
if (rustc.error) throw rustc.error;
if (rustc.status !== 0) {
  process.stderr.write(rustc.stderr);
  process.exit(rustc.status ?? 1);
}

const target = rustc.stdout.match(/^host: (.+)$/m)?.[1];
if (!target) throw new Error("Could not determine the Rust host target.");

run("cargo", ["build", "--release"], repoDir);

const executable = process.platform === "win32" ? "skiver.exe" : "skiver";
const sidecar = `skiver-${target}${process.platform === "win32" ? ".exe" : ""}`;
const destinationDir = join(uiDir, "src-tauri", "binaries");
mkdirSync(destinationDir, { recursive: true });
copyFileSync(
  join(repoDir, "target", "release", executable),
  join(destinationDir, sidecar),
);

console.log(`Prepared Tauri sidecar: ${sidecar}`);
