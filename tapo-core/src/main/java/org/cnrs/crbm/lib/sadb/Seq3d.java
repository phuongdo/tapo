package org.cnrs.crbm.lib.sadb;

import java.util.ArrayList;
import java.util.List;

/**
 * 
 * This class allows to combine 2, 3 or more Conformation alphabets into one.
 * This function helps reduce the long continuous identical alphabets which are
 * adjacent to each others. Hence, running TRUST and T-Recks are more efficient.
 * 
 * 
 * @author pdoviet
 *
 */
public class Seq3d {

	public static void main(String[] args) {
		String x1a50 = "xafaaaaaafaaglbvvvvsbvvsivgbaaafaaffaafaagxpgvvvvbvsvaxxhgbsaafaaffaafagglsbaaafaaffaafaaabggpsvvvbvvffgfaaxxaaaafaafaggxsavvvbaxbsfaghaafaaffaaxlvpsvbbvsaxpbaaafaafagapbevvpbpgbpebpbaxyyaafaafaaglppppvpvxlsgbaaafaafaagxsgvvvvaaaffaafagxagbaaafaaafaaafaaffaghx";
		String x1a4y = "xbxbvbvbvblbpfxaaafaafagahggsgbvbxvlhlsbaaghaafaaafagvggsgbvvxplvxfeaaffaafhghhxbxgpxsabvvxplhlsbaghhaaffaffggxgasgbvbhslvxfeaaffaafhafagpggpxsgbvbxplhlsbgahfaaffaffggbgxsabvvxslvxheaaffaaffaffaasgsxsabvbxplhlsbaahfaafaaffaapggsgbvvxslvxfeaaffaaffaffgpagpxsavvbxslhlvbaahfaaffaffggbggsabvvhslvxfegaffaaffafhgbpxsxsabvvxplhlsbaahfaafaaffaapgasgbvbhslvxhegaffaafaaffgpplpgsavvvxslhlsbgahhaafaaagaasgasgbvbhslbxsbaaffaafaaafapggvgsgbvvxslbpvpaafaaafgaaaaaxggvbbpx";
		String x1lxa = "xbpagpbbpagpbbppxsbvpplsbvsplvbvpggsbvpplvbvablvbvbexvbvxblvbvpplsbvippbgpgapxlpxsbvbvxblsbvpplsbvbpsvaagxlbvbvxblvbvpplsbvpplsbvxblvbvpplsbvpplvbvbplsbvpplvbvpplsbvsplsbvxxlvvvgbsvsslvbvbexxpbpgepxaaafaagxppaafaaffaafaaafgaalbbaaafaaafaaaggabaafaafaafagapagebsx";
		String xfake = "lxismmmmmbapxxxmmmmmm";
		Seq3d seq3d = new Seq3d(xfake);
		System.out.println(xfake);
		System.out.println(seq3d.getSeq());

	}

	public String getSeq() {
		StringBuilder builder = new StringBuilder();
		for (Conformation confor : confors) {
			builder.append(confor.getLetter());
		}
		return builder.toString();

	}

	List<Conformation> confors = new ArrayList<Conformation>();

	public List<Conformation> getConfors() {
		return confors;
	}

	public Seq3d(String strSA) {

		convertToConfor(strSA);

	}

	private void convertToConfor(String strSA) {

		int winsize = 2;

		String a_pattern = "";
		String f_pattern = "";
		String b_pattern = "";
		String v_pattern = "";
		for (int i = 0; i < winsize; i++) {
			a_pattern += "a";
			f_pattern += "f";
			b_pattern += "b";
			v_pattern += "v";
		}

		for (int i = 0; i < strSA.length() - winsize; i++) {

			String sub = strSA.substring(i, i + winsize);
			if (sub.equals(a_pattern)) {
				// o
				confors.add(new Conformation('o', i, i + winsize - 1));
				i = i + winsize - 1;
			} else if (sub.equals(v_pattern)) {
				// w
				confors.add(new Conformation('w', i, i + winsize - 1));
				i = i + winsize - 1;
			} else if (sub.equals(b_pattern)) {
				// m
				confors.add(new Conformation('m', i, i + winsize - 1));
				i = i + winsize - 1;
			} else if (sub.equals(f_pattern)) {

				// n
				confors.add(new Conformation('n', i, i + winsize - 1));
				i = i + winsize - 1;
			} else {

				confors.add(new Conformation(strSA.charAt(i), i, i));

			}

		}

	}
}
