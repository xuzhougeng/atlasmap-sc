import { cp, mkdir, readFile, rm, writeFile } from "node:fs/promises";
import { resolve } from "node:path";

const root = resolve(import.meta.dirname, "..");
const distDir = resolve(root, "dist");
const cfDir = resolve(root, "cloudflare");

function getBuildUpstreamOrigin() {
  const raw =
    process.env.UPSTREAM_ORIGIN ??
    process.env.UPSTREAM_HOST ??
    process.env.UPSTREAM;

  if (typeof raw !== "string" || raw.trim().length === 0) {
    throw new Error(
      "build:cloudflare requires UPSTREAM_ORIGIN (or UPSTREAM_HOST). " +
        "Example: UPSTREAM_ORIGIN=http://your.backend.server.com npm run build:cloudflare"
    );
  }

  const v = raw.trim();
  if (v.startsWith("http://") || v.startsWith("https://")) return v;
  return `http://${v}`;
}

async function main() {
  await mkdir(distDir, { recursive: true });

  // Cloudflare Pages expects `functions/` to be sibling of build output directory.
  // So we generate: frontend/functions/**  (same level as frontend/dist/**)
  const srcFunctions = resolve(cfDir, "functions");
  const dstFunctions = resolve(root, "functions");
  await rm(dstFunctions, { recursive: true, force: true });
  await cp(srcFunctions, dstFunctions, { recursive: true, force: true });

  // Inject upstream into functions/_proxy.js (no default baked in)
  const upstreamOrigin = getBuildUpstreamOrigin();
  const proxyPath = resolve(dstFunctions, "_proxy.js");
  const proxyContent = await readFile(proxyPath, "utf8");
  const replaced = proxyContent.replace(
    '"__BUILD_UPSTREAM_ORIGIN__"',
    JSON.stringify(upstreamOrigin)
  );
  if (replaced === proxyContent) {
    throw new Error(
      "Cloudflare proxy template placeholder not found; cannot inject upstream."
    );
  }
  await writeFile(proxyPath, replaced);

  // Copy Pages routing config into dist root (dist/_routes.json)
  await cp(resolve(cfDir, "_routes.json"), resolve(distDir, "_routes.json"), {
    force: true,
  });

  // Optional redirects: keep trailing slash behavior consistent
  try {
    await cp(resolve(cfDir, "_redirects"), resolve(distDir, "_redirects"), {
      force: true,
    });
  } catch {
    // ignore if not provided
  }

  // Small marker to make it obvious this dist is CF-ready
  await writeFile(resolve(distDir, ".cloudflare-pages"), "functions+routes\n");
}

await main();
