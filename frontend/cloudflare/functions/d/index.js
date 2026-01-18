import { proxyRequest } from "../_proxy.js";

export async function onRequest(context) {
  const url = new URL(context.request.url);
  // Normalize /d -> /d/
  if (url.pathname === "/d") {
    url.pathname = "/d/";
    return Response.redirect(url.toString(), 301);
  }

  return proxyRequest({
    request: context.request,
    context,
    // If you ever need "/d/*" -> upstream "/ad/*", enable:
    // transformPathname: (p) => p.replace(/^\/d(\/|$)/, "/ad$1"),
  });
}

