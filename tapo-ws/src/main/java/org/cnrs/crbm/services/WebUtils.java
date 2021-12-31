package org.cnrs.crbm.lib.services;

import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.util.HashMap;
import java.util.Map;

public class WebUtils {
	public static Map<String, String> getParameters(String queryString)
			throws UnsupportedEncodingException {
		Map<String, String> parameters = new HashMap<String, String>();
		String decoded = URLDecoder.decode(queryString, "UTF-8");
		String[] pares = decoded.split("&");
		for (String pare : pares) {
			String[] nameAndValue = pare.split("=");
			parameters.put(nameAndValue[0], nameAndValue[1]);
		}

		return parameters;
	}
}
