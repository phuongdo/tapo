package org.cnrs.crbm.lib.utils;

import org.xml.sax.SAXException;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.FileOutputStream;
import java.io.IOException;

//import org.fit.cssbox.css.CSSNorm;
//import org.fit.cssbox.css.DOMAnalyzer;
//import org.fit.cssbox.css.DOMAnalyzer.Origin;
//import org.fit.cssbox.demo.ImageRenderer;
//import org.fit.cssbox.io.DOMSource;
//import org.fit.cssbox.io.DefaultDOMSource;
//import org.fit.cssbox.io.DefaultDocumentSource;
//import org.fit.cssbox.io.DocumentSource;
//import org.fit.cssbox.layout.BrowserCanvas;
//import org.fit.cssbox.layout.BrowserConfig;
//import org.fit.cssbox.layout.Viewport;
//import org.fit.cssbox.render.SVGRenderer;
//import org.w3c.dom.Document;

public class ImageRenderPng {
	private boolean loadImages = true;
	private boolean loadBackgroundImages = false;

	public BufferedImage renderURL(String htmlContent) throws IOException,
			SAXException {

//		DocumentSource docSource = new ImageDocumentSource(htmlContent);
//		DOMSource parser = new DefaultDOMSource(docSource);
//		Document doc = parser.parse();
//		DOMAnalyzer da = new DOMAnalyzer(doc, docSource.getURL());
//		da.attributesToStyles();
//		da.addStyleSheet(null, CSSNorm.stdStyleSheet(),
//				DOMAnalyzer.Origin.AGENT);
//		da.addStyleSheet(null, CSSNorm.userStyleSheet(),
//				DOMAnalyzer.Origin.AGENT);
//		da.addStyleSheet(null, CSSNorm.formsStyleSheet(),
//				DOMAnalyzer.Origin.AGENT);
//		da.getStyleSheets();
//
//		BrowserCanvas contentCanvas = new BrowserCanvas(da.getRoot(), da,
//				docSource.getURL());
//		contentCanvas.getConfig().setLoadImages(this.loadImages);
//		contentCanvas.getConfig().setLoadBackgroundImages(
//				this.loadBackgroundImages);
//		contentCanvas.createLayout(new Dimension(1, 1));
//		docSource.close();
//
//		return contentCanvas.getImage();
		return null;
	}

	public static void main(String[] args) throws Exception {

		FileOutputStream os = new FileOutputStream("hello-world.png");
		ImageRenderPng r = new ImageRenderPng();
		BufferedImage bi = r.renderURL("http://www.google.fr");
		ImageIO.write(bi, "png", os);

		os.close();
	}

}
