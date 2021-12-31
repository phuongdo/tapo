package org.cnrs.crbm.lib.dssp;

import java.util.HashMap;
import java.util.Map;

/**
 * 
 * ASAView: Solvent Accessibility Graphics for proteins. Shandar Ahmad, M.
 * Michael Gromiha, Hamed Fawareh and Akinori Sarai BMC Bioinformatics (2004)
 * PUBLIC http://www.ncbi.nlm.nih.gov/pmc/articles/PMC420234/ WEB SERVER
 * http://www.abren.net/asaview/
 * ALY_X_ALY pattern
 * 
 * 
 * @author phuongdo
 * 
 */
public class RelativeASA {
	// /These values are (in Å2) 110.2 (Ala), 144.1 (Asp),140.4 (Cys), 174.7
	// (Glu), 200.7 (Phe), 78.7 (Gly), 181.9 (His), 185.0 (Ile), 205.7 (Lys),
	// 183.1 (Leu), 200.1 (Met),146.4 (Asn), 141.9 (Pro), 178.6 (Gln), 229.0
	// (Arg), 117.2 (Ser), 138.7 (Thr), 153.7 (Val), 240.5 (Trp), and 213.7
	// (Tyr) respectively.

	public static Map<String, Double> getAASurface() {
		Map<String, Double> map = new HashMap<String, Double>();

		map.put("A", 110.2);// Ala
		map.put("D", 144.1);// Asp
		map.put("C", 140.4);// Cys
		map.put("E", 174.7);// Glu
		map.put("F", 200.7);// Phe
		map.put("G", 78.7);// Gly
		map.put("H", 181.9);// His
		map.put("I", 185.0);// Ile
		map.put("K", 205.7); // Lys
		map.put("L", 183.1);// Leu
		map.put("M", 200.1); // Met
		map.put("N", 146.4);// Asn
		map.put("P", 141.9);// Pro
		map.put("Q", 178.6);// Gln
		map.put("R", 229.0);// Arg
		map.put("S", 117.2);// Ser
		map.put("T", 138.7);// Thr
		map.put("V", 153.7);// Val
		map.put("W", 240.5);// Trp
		map.put("Y", 213.7);// Tyr

		return map;

	}
}
