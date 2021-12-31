package apps;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.Row;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;
import java.util.List;

public class GeneratePdbImage {
	public static void main(String[] args) throws Exception {
		String input = "input/Repeat3DReport.csv";
		String outDir = "C:/Users/CRBM/Desktop/pdbimages";
		ProteinCSVReader csvReader = new ProteinCSVReader();
		List<Row> rows = csvReader.getData(input);
		StringBuilder builder = new StringBuilder();
		builder.append("<table>");
		builder.append("<tr><td>protein</td><td>Tag</td><td>3D Finder</td><tr/>");

		for (Row row : rows) {
			// System.out.println(strLine);
			String pdbCode = row.getPdbCode();
			String pdbChain = row.getPdbChain();
			System.out.println(pdbCode + "_" + pdbChain);

			// String imageUrl = "http://www.pdb.org/pdb/images/"
			// + pdbCode.toUpperCase() + "_asr_r_500.jpg";
			// String destinationFile = pdbCode.toUpperCase() + ".jpg";
			// saveImage(imageUrl, outDir + "/" + destinationFile);

			// build html file

			// tag builder
			String tag = "";
			if (row.getByEyeTot() == 1)// repeat
				tag = "<td valign='top' bgcolor='green' >1</td>";
			else if (row.getByEyeTot() == 2)// not repeat in 3D
				tag = "<td valign='top' bgcolor='red' >2</td>";
			else if (row.getByEyeTot() == 3)// not repeat in 3D
				tag = "<td valign='top' bgcolor='lightgray' >3</td>";

			else if (row.getByEyeTot() == 4)// not repeat in 3D
				tag = "<td valign='top' bgcolor='pink' >4</td>";
			else if (row.getByEyeTot() == 5)// not repeat in 3D
				tag = "<td valign='top' bgcolor='lightgreen' >5</td>";
			else if (row.getByEyeTot() == 6)// not repeat in 3D
				tag = "<td valign='top' bgcolor='lightgreen' >6</td>";
			else if (row.getByEyeTot() == 7)// not repeat in 3D
				tag = "<td valign='top' bgcolor='lightgreen' >7</td>";

			else
				// other
				tag = "<td>" + row.getByEyeTot() + "</td>";

			String conclusion = "";
			if (row.getConclusion() == 1)// repeat
				conclusion = "<td valign='top' bgcolor='green' >1</td>";
			else
				// not repeat in 3D
				conclusion = "<td valign='top' bgcolor='red' >2</td>";

			builder.append("<tr><td>" + pdbCode.toLowerCase() + "_" + pdbChain
					+ "</td>" + tag + conclusion + "<td><img src='"
					+ pdbCode.toUpperCase() + ".jpg' /></tr>");

		}

		builder.append("</table>");
		DataIO.writeToFile(builder.toString(), outDir + "/pdbview.html");
	}

	public static void saveImage(String imageUrl, String destinationFile)
			throws IOException {
		URL url = new URL(imageUrl);
		InputStream is = url.openStream();
		OutputStream os = new FileOutputStream(destinationFile);

		byte[] b = new byte[2048];
		int length;

		while ((length = is.read(b)) != -1) {
			os.write(b, 0, length);
		}

		is.close();
		os.close();
	}

}
