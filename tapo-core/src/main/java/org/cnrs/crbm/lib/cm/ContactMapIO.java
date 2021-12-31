package org.cnrs.crbm.lib.cm;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.SecondaryLoop;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import javax.imageio.ImageIO;

import org.cnrs.crbm.lib.conf.Dir;

public class ContactMapIO {

	public void saveCMtoFile(ContactMap cm) {

		String filename = cm.getPdbCode();
		String outputfile = Dir.CM_OUTPUT_DIR + "/" + filename + ".cm";
		FileOutputStream out;
		try {
			out = new FileOutputStream(outputfile);

			PrintStream p = new PrintStream(out);
			// System.out.println(cm.toFile());
			p.println(cm.toFile());
			p.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void saveMatrixtoPNG(int[][] cm, int nRows, int nCols) {
		try {
			String outputfile = Dir.CM_OUTPUT_DIR + "/pairsalign.png";

			int width = nCols, height = nRows;
			// TYPE_INT_ARGB specifies the image format: 8-bit RGBA packed
			// into integer pixels
			BufferedImage bi = new BufferedImage(width, height,
					BufferedImage.TYPE_INT_ARGB);
			Graphics2D ig2 = bi.createGraphics();

			// Font font = new Font("TimesRoman", Font.BOLD, 20);
			// ig2.setFont(font);
			ig2.setPaint(Color.black);

			for (int i = 0; i < nRows; i++) {

				for (int j = 0; j < nCols; j++) {

					if (cm[i][j] == 1) {
						ig2.drawLine(i, j, i, j);
					}
				}

			}

			ImageIO.write(bi, "PNG", new File(outputfile));

		} catch (IOException ie) {
			ie.printStackTrace();
		}
	}

	public void saveHistogramtoPNG(ContactMap cm) {
		try {

			int SIZE = 10;
			int BOTTON_MARGIN = 5;
			double[] histogram = cm.getHistogram();

			String outputfile = Dir.CM_OUTPUT_DIR + "/" + cm.getPdbCode()
					+ ".histo.png";

			int imgWidth = histogram.length, imgHeight = 10;

			for (int i = 0; i < histogram.length; i++) {
				if (imgHeight < histogram[i]) {
					imgHeight = (int) histogram[i];
				}
			}

			imgHeight = (imgHeight + BOTTON_MARGIN) * SIZE;
			imgWidth = imgWidth * SIZE;

			// TYPE_INT_ARGB specifies the image format: 8-bit RGBA packed
			// into integer pixels

			BufferedImage bi = new BufferedImage(imgWidth, imgHeight,
					BufferedImage.TYPE_INT_ARGB);
			Graphics2D ig2 = bi.createGraphics();

			// ig2.setFont(fonte);
			// ig2.drawString("bonjour", 50, 50);
			//
			ig2.setPaint(Color.gray);
			for (int i = 0; i < histogram.length; i++) {
				int maxHeight = (int) histogram[i];
				int x = i * SIZE;
				int y = (imgHeight / SIZE - BOTTON_MARGIN - maxHeight) * SIZE;
				ig2.fillRect(x, y, SIZE, maxHeight * SIZE);
			}
			// Font fonte = new Font("TimesRoman ", Font.BOLD, 10);
			for (int i = 0; i < histogram.length; i++) {
				if (cm.getSecondaryStructure().charAt(i) == 'H') {
					ig2.setColor(Color.red);
					ig2.fillRect(i * SIZE, imgHeight - BOTTON_MARGIN * SIZE,
							SIZE, BOTTON_MARGIN * SIZE);
				} else if (cm.getSecondaryStructure().charAt(i) == 'B') {
					ig2.setColor(Color.yellow);
					ig2.fillRect(i * SIZE, imgHeight - BOTTON_MARGIN * SIZE,
							SIZE, BOTTON_MARGIN * SIZE);
				}
			}

			ImageIO.write(bi, "PNG", new File(outputfile));

		} catch (IOException ie) {
			ie.printStackTrace();
		}
	}

	public void saveCMtoPNG(ContactMap cm) {
		try {
			String filename = cm.getPdbCode();
			String outputfile = Dir.CM_OUTPUT_DIR + "/" + filename + ".png";

			int width = cm.getNoRes(), height = cm.getNoRes();
			// TYPE_INT_ARGB specifies the image format: 8-bit RGBA packed
			// into integer pixels
			BufferedImage bi = new BufferedImage(width, height,
					BufferedImage.TYPE_INT_ARGB);
			Graphics2D ig2 = bi.createGraphics();

			// Font font = new Font("TimesRoman", Font.BOLD, 20);
			// ig2.setFont(font);
			ig2.setPaint(Color.black);
			for (Pair pair : cm.getCmaps()) {

				int x = pair.getFirst();
				int y = pair.getSecond();
				ig2.drawLine(x, y, x, y);
				ig2.drawLine(y, x, y, x);
			}

			ImageIO.write(bi, "PNG", new File(outputfile));

		} catch (IOException ie) {
			ie.printStackTrace();
		}
	}

	public void readCMfromFile() {
	}

}
