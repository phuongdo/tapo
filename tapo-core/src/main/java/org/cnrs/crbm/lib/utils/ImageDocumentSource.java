//package org.cnrs.crbm.lib.utils;
//
//import java.io.ByteArrayInputStream;
//import java.io.IOException;
//import java.io.InputStream;
//import java.net.URL;
//import java.net.URLConnection;
//
//
////import org.fit.cssbox.io.DocumentSource;
////import org.fit.net.DataURLHandler;
//
//public class ImageDocumentSource extends DocumentSource {
//	private static String USER_AGENT = "Mozilla/5.0 (compatible; BoxBrowserTest/4.x; Linux) CSSBox/4.x (like Gecko)";
//	private URLConnection con;
//	private InputStream is;
//	private String htmlcontent = "";
//
//	public ImageDocumentSource(String htmlContent) throws IOException {
//		super(null, null);
//		this.con = null;
//		this.is = null;
//		this.htmlcontent = htmlContent;
//	}
//
//	public ImageDocumentSource(URL url) throws IOException {
//		super(url);
//		this.con = createConnection(url);
//		this.is = null;
//	}
//
//	// public ImageDocumentSource(String urlstring) throws IOException {
//	// super(null, urlstring);
//	// URL url = DataURLHandler.createURL(null, urlstring);
//	// this.con = createConnection(url);
//	// this.is = null;
//	// }
//
//	public ImageDocumentSource(URL base, String urlstring) throws IOException {
//		super(base, urlstring);
//		URL url = DataURLHandler.createURL(base, urlstring);
//		this.con = createConnection(url);
//		this.is = null;
//	}
//
//	protected URLConnection createConnection(URL url) throws IOException {
//		URLConnection con = url.openConnection();
//		con.setRequestProperty("User-Agent", USER_AGENT);
//		return con;
//	}
//
//	public URL getURL() {
//		// return this.con.getURL();
//
//		return null;
//
//	}
//
//	public InputStream getInputStream() throws IOException {
//		if (this.is == null) {
//			this.is = new ByteArrayInputStream(this.htmlcontent.getBytes());
//		}
//		return this.is;
//	}
//
//	public String getContentType() {
//		// return this.con.getHeaderField("Content-Type");
//		return "none";
//	}
//
//	public static String getUserAgent() {
//		return USER_AGENT;
//	}
//
//	public static void setUserAgent(String userAgent) {
//		USER_AGENT = userAgent;
//	}
//
//	public void close() throws IOException {
//		if (this.is != null) {
//			this.is.close();
//		}
//	}
//}
