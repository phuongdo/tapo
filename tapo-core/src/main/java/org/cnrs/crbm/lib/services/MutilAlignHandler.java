package org.cnrs.crbm.lib.services;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Map;

import javax.imageio.ImageIO;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.handler.AbstractHandler;

/**
 * http://localhost:8888/mutilalign/?code=1LXA&chain=A&align=30-39;42-51;54-63
 * 
 * @author phuongdo
 * 
 */
public class MutilAlignHandler extends AbstractHandler {

	public void handle(String target, Request baseRequest,
			HttpServletRequest request, HttpServletResponse response)
			throws IOException, ServletException {
		try {
			response.setContentType("text/html;charset=utf-8");
			// response.setContentType("image/png");
			response.setStatus(HttpServletResponse.SC_OK);
			baseRequest.setHandled(true);

			String queryString = request.getQueryString();
			Map<String, String> parameters = WebUtils
					.getParameters(queryString);

			String code = parameters.get("code").toLowerCase();
			String chain = parameters.get("chain").toUpperCase();
			String align = parameters.get("align");

			MutilAlign mutilAlig = new MutilAlign(code, chain, align);
			response.getWriter().write(mutilAlig.getOutputHtml());

			// BufferedImage bufferedImage = mutilAlig.getImage();
			// OutputStream out = response.getOutputStream();
			// ImageIO.write(bufferedImage, "png", out);

			// out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
