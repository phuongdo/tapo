package org.cnrs.crbm.lib.io;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.zip.GZIPInputStream;

public class GzipIO {

	public static void extract(String gzipFile, String outputFile) {
		byte[] buffer = new byte[1024];
		try {

			GZIPInputStream gzis = new GZIPInputStream(new FileInputStream(
					gzipFile));

			FileOutputStream out = new FileOutputStream(outputFile);

			int len;
			while ((len = gzis.read(buffer)) > 0) {
				out.write(buffer, 0, len);
			}
			gzis.close();
			out.close();
			// System.out.println("Done");

		} catch (IOException ex) {
			ex.printStackTrace();
		}

	}

}
