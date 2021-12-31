package demo;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class NaidaProblem {

	public static Map<String, ArrayList<String>> getParentHeader() {
		Map<String, ArrayList<String>> map = new HashMap<String, ArrayList<String>>();
		try {
			// Open the file that is the first
			// command line parameter
			FileInputStream fstream = new FileInputStream("input/nadia");
			// Get the object of DataInputStream
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String line;

			// Read File Line By Line

			while ((line = br.readLine()) != null) {

				// extract header

				String header = line.substring(0, 2);
				String content = line.substring(3, line.length());

				if (map.containsKey(header)) {
					map.get(header).add(content);

				} else {

					ArrayList<String> list = new ArrayList<String>();
					list.add(content);
					map.put(header, list);

				}

			}

			// Close the input stream
			in.close();
		} catch (Exception e) {// Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}

		return map;

	}

	public static void main(String[] args) {

		Map<String, ArrayList<String>> map = getParentHeader();

		ArrayList<String> subContent = map.get("CC");

		Map<Integer, String> headers = new HashMap<Integer, String>();
		Map<Integer, String> content = new HashMap<Integer, String>();
		int hindex = 0;
		for (String str : subContent) {

			// System.out.println(str);
			str = str.trim();
			if (str.equals(""))
				continue;

			if (str.startsWith("-!-")) {
				hindex++;
				headers.put(hindex, str.substring(0, 8));
				content.put(hindex, str.substring(8, str.length() - 1));

			} else {
				if (content.containsKey(hindex))
					content.put(hindex, content.get(hindex) + str);

			}

			// extract sub-content header;
		}

		System.out.println(headers.get(8));
		System.err.println(content.get(8));
	}
}
