function getUpstreamOrigin(context) {
  // Prefer Pages env var if configured; fall back to default.
  // In Cloudflare Pages: Settings -> Environment variables.
  const fromEnv =
    context?.env?.UPSTREAM_ORIGIN ||
    context?.env?.UPSTREAM_HOST ||
    context?.env?.UPSTREAM;

  if (typeof fromEnv === "string" && fromEnv.trim().length > 0) {
    const v = fromEnv.trim();
    // Allow either "wanglab.sippe.ac.cn" or "http://wanglab.sippe.ac.cn"
    if (v.startsWith("http://") || v.startsWith("https://")) return v;
    return `http://${v}`;
  }

  // This value is injected at build time by `npm run build:cloudflare`.
  // It intentionally has NO baked-in default.
  const injected = "__BUILD_UPSTREAM_ORIGIN__";
  if (injected && injected !== "__BUILD_UPSTREAM_ORIGIN__") return injected;

  throw new Error(
    "Missing upstream origin. Set UPSTREAM_ORIGIN/UPSTREAM_HOST at build time " +
      "or configure UPSTREAM_ORIGIN/UPSTREAM_HOST in Cloudflare Pages environment variables."
  );
}

export async function proxyRequest({ request, context, transformPathname }) {
  const publicUrl = new URL(request.url);

  const upstreamBase = new URL(getUpstreamOrigin(context));
  const upstreamUrl = new URL(publicUrl.toString());
  upstreamUrl.protocol = upstreamBase.protocol;
  upstreamUrl.host = upstreamBase.host;

  if (typeof transformPathname === "function") {
    upstreamUrl.pathname = transformPathname(upstreamUrl.pathname);
  }

  const init = {
    method: request.method,
    headers: new Headers(request.headers),
    body: request.body,
    redirect: "manual",
  };

  // Let upstream know the user-facing scheme + host.
  init.headers.set("X-Forwarded-Proto", publicUrl.protocol.replace(":", ""));
  init.headers.set("X-Forwarded-Host", publicUrl.host);

  const resp = await fetch(upstreamUrl.toString(), init);
  return rewriteLocation(resp, upstreamUrl, publicUrl);
}

function rewriteLocation(resp, upstreamUrl, publicUrl) {
  const r = new Response(resp.body, resp);
  const loc = r.headers.get("Location");
  if (!loc) return r;

  try {
    const locUrl = new URL(loc, upstreamUrl);
    // If upstream redirects back to itself, rewrite to same-origin.
    if (locUrl.host === upstreamUrl.host) {
      const newLoc = new URL(publicUrl.origin);
      newLoc.pathname = locUrl.pathname;
      newLoc.search = locUrl.search;
      r.headers.set("Location", newLoc.toString());
    }
  } catch {
    // ignore parse errors
  }

  return r;
}

