package org.cnrs.crbm.lib.services;

import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Map;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.randompdb.RandomPDB;
import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.handler.AbstractHandler;

/**
 * http://localhost:8888/fakestruct/
 * 
 * @author phuongdo
 * 
 */
public class FakeStructureHandler extends AbstractHandler {

	public void handle(String target, Request baseRequest,
			HttpServletRequest request, HttpServletResponse response)
			throws IOException, ServletException {
		try {
			response.setContentType("text/html;charset=utf-8");
			response.setStatus(HttpServletResponse.SC_OK);
			baseRequest.setHandled(true);

			String queryString = request.getQueryString();

			RandomPDB random = new RandomPDB();
			response.getWriter().println(random.build());

		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
