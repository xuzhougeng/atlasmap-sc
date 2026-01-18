import { proxyRequest } from "../_proxy.js";

export async function onRequest(context) {
  return proxyRequest({
    request: context.request,
    context,
    // If you ever need "/d/*" -> upstream "/ad/*", enable:
    // transformPathname: (p) => p.replace(/^\/d(\/|$)/, "/ad$1"),
  });
}

