package org.cnrs.crbm.lib.annotation;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.trsfinder.Features;

public class ProteinAnnotate {

	/**
	 * sequence repeats are in one secondary structure, usually in long Helix as
	 * static, Helix length > 40 residues or > 1/2 length of this protein is
	 * considered as protein with long Helix
	 * 
	 * @return
	 */

	private final static String MESG_LONG_HELIX = "LH";// "1.Long Helix";
	private final static String MESG_MISSING_ATOMS = "MA";// "2.Missing Atoms";
	private final static String MESG_GLOBULAR = "GL";// "3.Globular";

	public static String anotation(RepeatFinder finder, int start, int end) {

		ProteinAnnotate proAnotate = new ProteinAnnotate();

		StringBuilder builder = new StringBuilder();

		if (proAnotate.annotateProteinWithLongHelix(finder.getStrSS(), start,
				end)) {

			builder.append(MESG_LONG_HELIX + " ");
		}

		if (proAnotate.annotateDispersedProtein(finder.getAtoms(), start, end)) {
			builder.append(MESG_MISSING_ATOMS + " ");
		}

		if (proAnotate
				.annotateGlobularProtein(finder.getCmHistos(), start, end)) {
			builder.append(MESG_GLOBULAR + " ");
		}

		return builder.toString();

	}

	public static String anotation(Features feature, int start, int end) {

		ProteinAnnotate proAnotate = new ProteinAnnotate();

		StringBuilder builder = new StringBuilder();

		if (proAnotate.annotateProteinWithLongHelix(feature.getStrSS(), start,
				end)) {

			builder.append(MESG_LONG_HELIX + " ");
		}

		if (proAnotate.annotateDispersedProtein(feature.getAtoms(), start, end)) {
			builder.append(MESG_MISSING_ATOMS + " ");
		}

		if (proAnotate.annotateGlobularProtein(feature.getCmHistos(), start,
				end)) {
			builder.append(MESG_GLOBULAR + " ");
		}

		return builder.toString();

	}

	public boolean annotateProteinWithLongHelix(String secondaryStr,
			int startRepeat, int endRepeat) {

		boolean isLongHelix = false;
		int start = 0;
		int end = 0;
		int sizeRepeat = endRepeat - startRepeat + 1;
		for (int i = startRepeat; i < endRepeat; i++) {
			char ch = secondaryStr.charAt(i);
			if (ch == 'B') {
				start = i;
				while (i < secondaryStr.length()
						&& secondaryStr.charAt(i) == 'B') {
					i++;
				}
				i = end = i - 1;

			} else if (ch == 'H') {

				start = i;
				while (i < secondaryStr.length()
						&& secondaryStr.charAt(i) == 'H') {
					i++;
				}
				i = end = i - 1;

				int size = end - start + 1;

				if (size > 0.8 * sizeRepeat) {
					isLongHelix = true;
					break;
				}

			}

		}

		return isLongHelix;

	}

	/**
	 * deal with Rossmann fold like globular protein, too many other structures
	 * surrounded
	 * 
	 * @return
	 */
	public boolean annotateGlobularProtein(double[] cmHistos, int start, int end) {

		boolean isGlobular = false;
		int size = end - start + 1;
		// double threshold = (1 / 3) * cmHistos.length;
		double threshold = size;
		int count = 0;
		for (int i = start; i < end; i++) {
			if (cmHistos[i] > threshold)
				count++;
		}
		// System.out.print(count + "\t" + cmHistos.length + "\n");

		// System.out.println(count + ":" + 0.8 * size);
		if (count > 0.8 * size)
			isGlobular = true;

		return isGlobular;
	}

	private final static int THRESHOLD_DISPERSED = 0;

	/**
	 * annotate protein which is missing atom
	 * 
	 * @param atoms
	 * @return
	 */
	public boolean annotateDispersedProtein(Atom[] atoms, int start, int end) {

		// System.out.println(atoms.length + ":" + start + ":" + end);
		boolean isDispersed = false;
		int nrDispersed = 0;
		for (int i = start; i < end - 1; i++) {
			int seqCurr = atoms[i].getGroup().getResidueNumber().getSeqNum();
			int seqNext = atoms[i + 1].getGroup().getResidueNumber()
					.getSeqNum();
			if (seqCurr + 1 != seqNext)
				nrDispersed++;
		}

		// System.out.println(nrDispersed);
		if (nrDispersed > THRESHOLD_DISPERSED)
			isDispersed = true;

		return isDispersed;

	}
}
