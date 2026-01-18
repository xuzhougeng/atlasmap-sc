import { proxyRequest } from "../_proxy.js";

export async function onRequest(context) {
  const url = new URL(context.request.url);
  // Normalize /api -> /api/
  if (url.pathname === "/api") {
    url.pathname = "/api/";
    return Response.redirect(url.toString(), 301);
  }

  return proxyRequest({ request: context.request, context });
}

