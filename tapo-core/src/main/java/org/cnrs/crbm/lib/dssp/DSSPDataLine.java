package org.cnrs.crbm.lib.dssp;

import java.util.Locale;

public class DSSPDataLine {
	public int position;
	public int seqPosition;
	public char chainLetter;
	public char residueLetter;
	public char structure;
	public int bp1;
	public int bp2;
	public int acc;
	public int nhoPartner1;
	public double nhoEnergy1;
	public int ohnPartner1;
	public double ohnEnergy1;
	public int nhoPartner2;
	public double nhoEnergy2;
	public int ohnPartner2;
	public double ohnEnergy2;
	public double tco;
	public double kappa;
	public double alpha;
	public double phi;
	public double psi;
	public double xCa;
	public double yCa;
	public double zCa;

	public DSSPDataLine(String paramString, boolean paramBoolean)
			throws NumberFormatException {
		int[] arrayOfInt1 = a;
		int[] arrayOfInt2 = b;
		if (paramBoolean) {
			arrayOfInt1 = c;
			arrayOfInt2 = d;
		}
		this.position = Integer.valueOf(
				paramString.substring(arrayOfInt1[0] - 1, arrayOfInt2[0] - 1)
						.trim()).intValue();

		this.seqPosition = Integer.valueOf(
				paramString.substring(arrayOfInt1[1] - 1, arrayOfInt2[1] - 1)
						.trim()).intValue();

		this.chainLetter = paramString.charAt(arrayOfInt1[2] - 1);

		this.residueLetter = paramString.charAt(arrayOfInt1[3] - 1);
		if (Character.isLowerCase(this.residueLetter)) {
			this.residueLetter = 'C';
		}
		try {
			this.structure = paramString.charAt(arrayOfInt1[4] - 1);
			this.bp1 = Integer.valueOf(
					paramString.substring(arrayOfInt1[5] - 1,
							arrayOfInt2[5] - 1).trim()).intValue();

			this.bp2 = Integer.valueOf(
					paramString.substring(arrayOfInt1[6] - 1,
							arrayOfInt2[6] - 1).trim()).intValue();

			this.acc = Integer.valueOf(
					paramString.substring(arrayOfInt1[7] - 1,
							arrayOfInt2[7] - 1).trim()).intValue();

			this.nhoPartner1 = Integer.valueOf(
					paramString.substring(arrayOfInt1[8] - 1,
							arrayOfInt2[8] - 1).trim()).intValue();

			this.nhoEnergy1 = Double.valueOf(
					paramString.substring(arrayOfInt1[9] - 1,
							arrayOfInt2[9] - 1).trim()).doubleValue();

			this.ohnPartner1 = Integer.valueOf(
					paramString.substring(arrayOfInt1[10] - 1,
							arrayOfInt2[10] - 1).trim()).intValue();

			this.ohnEnergy1 = Double.valueOf(
					paramString.substring(arrayOfInt1[11] - 1,
							arrayOfInt2[11] - 1).trim()).doubleValue();

			this.nhoPartner2 = Integer.valueOf(
					paramString.substring(arrayOfInt1[12] - 1,
							arrayOfInt2[12] - 1).trim()).intValue();

			this.nhoEnergy2 = Double.valueOf(
					paramString.substring(arrayOfInt1[13] - 1,
							arrayOfInt2[13] - 1).trim()).doubleValue();

			this.ohnPartner2 = Integer.valueOf(
					paramString.substring(arrayOfInt1[14] - 1,
							arrayOfInt2[14] - 1).trim()).intValue();

			this.ohnEnergy2 = Double.valueOf(
					paramString.substring(arrayOfInt1[15] - 1,
							arrayOfInt2[15] - 1).trim()).doubleValue();

			this.tco = Double.valueOf(
					paramString.substring(arrayOfInt1[16] - 1,
							arrayOfInt2[16] - 1).trim()).doubleValue();
		
		this.kappa = Double.valueOf(
				paramString.substring(arrayOfInt1[17] - 1, arrayOfInt2[17] - 1)
						.trim()).doubleValue();

		this.alpha = Double.valueOf(
				paramString.substring(arrayOfInt1[18] - 1, arrayOfInt2[18] - 1)
						.trim()).doubleValue();

		this.phi = Double.valueOf(
				paramString.substring(arrayOfInt1[19] - 1, arrayOfInt2[19] - 1)
						.trim()).doubleValue();

		this.psi = Double.valueOf(
				paramString.substring(arrayOfInt1[20] - 1, arrayOfInt2[20] - 1)
						.trim()).doubleValue();

		this.xCa = Double.valueOf(
				paramString.substring(arrayOfInt1[21] - 1, arrayOfInt2[21] - 1)
						.trim()).doubleValue();

		this.yCa = Double.valueOf(
				paramString.substring(arrayOfInt1[22] - 1, arrayOfInt2[22] - 1)
						.trim()).doubleValue();

		this.zCa = Double.valueOf(
				paramString.substring(arrayOfInt1[23] - 1, arrayOfInt2[23] - 1)
						.trim()).doubleValue();

		} catch (Exception ex) {

			System.out.println(paramString);
			// ex.printStackTrace();

		}
	}

	public String toString() {
		return String.format(
				Locale.ENGLISH,
				"%5d %4d %c %c %c %8.3f %8.3f %8.3f %8.3f %8.3f",
				new Object[] { Integer.valueOf(this.position),
						Integer.valueOf(this.seqPosition),
						Character.valueOf(this.chainLetter),
						Character.valueOf(this.residueLetter),
						Character.valueOf(this.structure),
						Double.valueOf(this.tco), Double.valueOf(this.alpha),
						Double.valueOf(this.kappa), Double.valueOf(this.phi),
						Double.valueOf(this.psi) });
	}

	private static final int[] a = { 1, 7, 12, 14, 17, 27, 31, 35, 40, 47, 52,
			58, 63, 69, 74, 80, 86, 92, 98, 104, 110, 116, 123, 130 };
	private static final int[] b = { 6, 11, 13, 15, 18, 30, 34, 39, 46, 51, 57,
			62, 68, 73, 79, 84, 92, 98, 104, 110, 116, 123, 130, 137 };
	private static final int[] c = { 1, 7, 12, 14, 17, 27, 31, 35, 40, 45, 49,
			54, 58, 63, 67, 72, 76, 84, 90, 96, 102, 108, 115, 122 };
	private static final int[] d = { 6, 11, 13, 15, 18, 30, 34, 39, 44, 49, 53,
			58, 62, 67, 71, 76, 84, 90, 96, 102, 108, 115, 120, 129 };
}