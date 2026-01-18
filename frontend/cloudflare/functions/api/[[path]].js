import { proxyRequest } from "../_proxy.js";

export async function onRequest(context) {
  return proxyRequest({ request: context.request, context });
}

