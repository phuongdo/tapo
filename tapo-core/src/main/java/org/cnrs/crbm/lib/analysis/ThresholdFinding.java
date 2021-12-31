package org.cnrs.crbm.lib.analysis;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.RowRepeatDB;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.sadb.Seq3d;
import org.cnrs.crbm.lib.seqalign.TREksWapper;
import org.cnrs.crbm.lib.trsfinder.Features;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class ThresholdFinding {
	DecimalFormat df = new DecimalFormat("0.00");

	public static void main(String[] args) throws Exception {

        // System.out.println("2-95".split("-")[0]);

        new ThresholdFinding().findQADistribution();
        //new ThresholdFinding().run();

	}


    void findQADistribution() throws Exception {

        ProteinCSVReader csvReader = new ProteinCSVReader();

        List<RowRepeatDB> rows = csvReader.getRepeatDB("data/RDB-dataset.tab");
        PrintWriter writer = new PrintWriter("data/threshold_QAx.train");
        writer.write("pdb\tclass\tavgL\tscore1\tscore2\tQA\n");
        // "UTF-8");
        for (RowRepeatDB row : rows) {

            //

            if (row.getAnnlevel().equals("Detailed") && row.getUnits().length() > 0 && row.getUnits().split(";").length < 10) {
                // System.out.println(row.getEntry());


                try {
                    String[] units = row.getUnits().split(";");
                    String pdbCode = row.getPdbCode();
                    String pdbChain = row.getPdbChain();
                    Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
                    Atom[] atoms = PdbTools.getAtomCAArray(structure
                            .getChainByPDB(pdbChain));


                    RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);

                    StringBuffer buffer = new StringBuffer();
                    double avgL = 0.0;
                    for (String unit : units) {
                        int starti = Integer.parseInt(unit.split("-")[0]);
                        int endi = Integer.parseInt(unit.split("-")[1]);

                        int start_pos = this.getPosition(atoms,
                                starti);
                        int end_pos = this.getPosition(atoms, endi);
                        avgL += end_pos - start_pos + 1;

                        buffer.append(start_pos + "-" + end_pos + ";");

                    }

                    avgL = avgL / units.length;
                    String frags = buffer.toString();
                    frags = frags.substring(0, frags.length() - 1);


                    // convert to
                    MutilAlign mutilAlign = new MutilAlign(repeatFinder.getFeatures(), frags);
                    writer.write(row.getEntry() + "\t" + row.getStrclass() + "\t" + NumberFormatUtils.format(avgL) + "\t" + NumberFormatUtils.format(mutilAlign.getScore1()) + "\t" + NumberFormatUtils.format(mutilAlign.getScore2()) + "\t" + NumberFormatUtils.format(mutilAlign.getQAScore()) + "\n");
                    System.out.println(row.getEntry() + "\t" + NumberFormatUtils.format(mutilAlign.getScore1()) + "\t" + NumberFormatUtils.format(mutilAlign.getScore2()) + "\t" + NumberFormatUtils.format(mutilAlign.getQAScore()));
                } catch (Exception ex) {
                    ex.printStackTrace();
                }

            }
        }
        writer.close();

    }



	void run() throws Exception {
		ProteinCSVReader csvReader = new ProteinCSVReader();

		List<RowRepeatDB> rows = csvReader.getRepeatDB("data/RDB-dataset.tab");
        PrintWriter writer = new PrintWriter("data/threshold.train");
        writer.write("note\taleng\tRMSD\tgapRes\tleng1\tleng2\tminL\tmaxL\taL\taS\ttmScore\tcolor\n");
        // "UTF-8");
        for (RowRepeatDB row : rows) {

            if (row.getAnnlevel().equals("Detailed") && row.getUnits().length() > 0) {

                //System.out.println(row.getPdbCode() + "_" + row.getPdbChain());
                // processTrust(row);
                // processTReks(row);
                process(row, writer);
                // break;
            }
        }
        writer.close();
    }

	void processTReks(RowRepeatDB row) throws StructureException {

		String pdbCode = row.getPdbCode();
		String pdbChain = row.getPdbChain();
		String[] units = row.getUnits().split(";");
		RepeatFinder finder = new RepeatFinder(pdbCode, pdbChain);
		try {
			finder.findByOneMethod();

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// Features features = finder.getFeatures();

		// Seq3d seq3d = new Seq3d(features.getStrSeqAp());

		// int gapo = -5;
		// int gape = -1;
		// Trust trust = new Trust("conf/BLOSUM_COMBINED", pdbCode,
		// seq3d.getSeq()
		// .toUpperCase(), gapo, gape);

		// Sequence3D sequence3D = new Sequence3D();
		// // String pattern = sequence3D.getSAV3(features.getStrSeqSADB());
		// String pattern = sequence3D.getSAV3(features.getStrSeqAp());
		// Seq3d seq3d = new Seq3d(pattern);
		// pattern = seq3d.getSeq();
		// TREksWapper TReks = new TREksWapper("", pattern);
		// // get treks repeats
		// List<nptr.Repeat> treks_repeats = TReks.getRepeats();
		// try {
		//
		// boolean isTRs = true;
		//
		// if (treks_repeats != null) {
		// int repeatsize = treks_repeats.size();
		//
		// if (repeatsize > 1) {
		// System.out.println(pdbCode + " : yes : " + repeatsize
		// + "  units:" + units.length);
		// } else {
		// System.out.println(pdbCode + " : no : " + repeatsize
		// + "  units:" + units.length);
		// isTRs = false;
		// }
		// } else {
		// System.out.println(pdbCode + " : no : ");
		//
		// isTRs = false;
		// }
		//
		// } catch (Exception ex) {
		// ex.printStackTrace();
		// }

	}

	void processTrust(RowRepeatDB row) throws StructureException {

		String pdbCode = row.getPdbCode();
		String pdbChain = row.getPdbChain();
		String[] units = row.getUnits().split(";");
		RepeatFinder finder = new RepeatFinder(pdbCode, pdbChain);

		Features features = finder.getFeatures();

		Seq3d seq3d = new Seq3d(features.getStrSeqAp());

		int gapo = -5;
		int gape = -1;
		// Trust trust = new Trust("conf/BLOSUM_COMBINED", pdbCode,
		// seq3d.getSeq()
		// .toUpperCase(), gapo, gape);

		TREksWapper TReks = new TREksWapper("", features.getStrSeqAp());
		// get treks repeats
		List<nptr.Repeat> treks_repeats = TReks.getRepeats();
		try {

			boolean isTRs = true;

			if (treks_repeats != null) {
				int repeatsize = treks_repeats.size();

				if (repeatsize > 1) {
					System.out.println(pdbCode + " : yes : " + repeatsize
							+ "  units:" + units.length);
				} else {
					System.out.println(pdbCode + " : no : " + repeatsize
							+ "  units:" + units.length);
					isTRs = false;
				}
			} else {
				System.out.println(pdbCode + " : no : ");

				isTRs = false;
			}

			// if (!isTRs) {
			//
			// // debug mode;
			// System.out.println(">>debug ici: ");
			// System.out.println(">>" + seq3d);
			//
			// for (int i = 0; i < units.length - 1; i++) {
			//
			// int starti = Integer.parseInt(units[i].split("-")[0]);
			// int endi = Integer.parseInt(units[i].split("-")[1]);
			//
			// int start_pos = this.getPosition(features.getAtoms(),
			// starti);
			// int end_pos = this.getPosition(features.getAtoms(), endi);
			// String alp = features.getStrSeqAp().substring(start_pos,
			// end_pos);
			//
			// Seq3d seq3dUnit = new Seq3d(alp);
			// System.out.println(">> " + seq3dUnit.getSeq());
			//
			// }
			//
			// }
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}


    public Set<String> tripeUnitsProcess(String[] units, Atom[] atoms, PrintWriter writer, RowRepeatDB row) {

        Set<String> hashSet = new HashSet<String>();

        MutilAlign mutilAlign = new MutilAlign();
        if (units.length >= 3) {

            int winSize = 3;
            try {


                for (int i = 0; i < units.length - winSize; i++) {

                    if (units[i].startsWith("-"))
                        continue;

                    for (int j = i; j < i + winSize; j++) {

                        int startj = Integer.parseInt(units[j].split("-")[0]);
                        int endj = Integer.parseInt(units[j].split("-")[1]);
                        Atom[] atomj = this.getAtoms(atoms, startj, endj);

                        for (int k = j + 1; k < i + winSize; k++) {

                            //System.out.println(units[i]);
                            int starti = Integer.parseInt(units[k].split("-")[0]);
                            int endi = Integer.parseInt(units[k].split("-")[1]);
                            Atom[] atomi = this.getAtoms(atoms, starti, endi);


                            AFPChain afpChain = null;
                            afpChain = mutilAlign.pairAlign(atomi, atomj);


                            double rmsd = afpChain.getChainRmsd();
                            int alnLeng = afpChain.getOptLength();
                            int gapLeng = afpChain.getGapLen();
                            int leng1 = afpChain.getCa1Length();
                            int leng2 = afpChain.getCa2Length();
                            int minL = Math.min(leng1, leng2);
                            int maxL = Math.max(leng1, leng2);
                            double aL = (double) alnLeng / maxL;
                            double aS = (double) alnLeng / minL;
                            double tmScore = afpChain.getTMScore();


                            if (aL < 0.9 || aS < 0.9 && tmScore > 0)
                                continue;
//                            writer.write(row.getEntry() + "(" + units[j] + " vs " + units[k] + ")" + "\t" + alnLeng + "\t" + df.format(rmsd) + "\t"
//                                    + gapLeng + "\t" + leng1 + "\t" + leng2 + "\t" + minL + "\t" + maxL + "\t" + df.format(aL) + "\t" + df.format(aS) + "\t"
//                                    + df.format(tmScore) + "\t" + getClassOfScore(gapLeng) + "\n");

                            hashSet.add(row.getEntry() + "(" + units[j] + " vs " + units[k] + ")" + "\t" + alnLeng + "\t" + df.format(rmsd) + "\t"
                                    + gapLeng + "\t" + leng1 + "\t" + leng2 + "\t" + minL + "\t" + maxL + "\t" + df.format(aL) + "\t" + df.format(aS) + "\t"
                                    + df.format(tmScore) + "\t" + getClassOfScore(gapLeng));
                        }
                    }
                }

            } catch (StructureException e) {
                //e.printStackTrace();
            }


        } else if (units.length == 2)

        {

            try {
                int starti = Integer.parseInt(units[0].split("-")[0]);
                int endi = Integer.parseInt(units[0].split("-")[1]);
                Atom[] atomi = this.getAtoms(atoms, starti, endi);
                int startj = Integer.parseInt(units[1].split("-")[0]);
                int endj = Integer.parseInt(units[1].split("-")[1]);

                Atom[] atomj = this.getAtoms(atoms, startj, endj);
                AFPChain afpChain = null;
                afpChain = mutilAlign.pairAlign(atomi, atomj);
                double rmsd = afpChain.getChainRmsd();
                int alnLeng = afpChain.getOptLength();
                int gapLeng = afpChain.getGapLen();
                int leng1 = afpChain.getCa1Length();
                int leng2 = afpChain.getCa2Length();
                int minL = Math.min(leng1, leng2);
                int maxL = Math.max(leng1, leng2);
                double aL = (double) alnLeng / maxL;
                double aS = (double) alnLeng / minL;

                double tmScore = afpChain.getTMScore();


//                writer.write(row.getEntry() + "(" + units[0] + " vs " + units[1] + ")" + "\t" + alnLeng + "\t" + df.format(rmsd) + "\t"
//                        + gapLeng + "\t" + leng1 + "\t" + leng2 + "\t" + minL + "\t" + maxL + "\t" + aL + "\t" + aS + "\t"
//                        + df.format(tmScore) + "\t" + getClassOfScore(gapLeng) + "\n");

                hashSet.add(row.getEntry() + "(" + units[0] + " vs " + units[1] + ")" + "\t" + alnLeng + "\t" + df.format(rmsd) + "\t"
                        + gapLeng + "\t" + leng1 + "\t" + leng2 + "\t" + minL + "\t" + maxL + "\t" + df.format(aL) + "\t" + df.format(aS) + "\t"
                        + df.format(tmScore) + "\t" + getClassOfScore(gapLeng));

            } catch (StructureException e) {
                e.printStackTrace();
            }
        }

        return hashSet;

    }

    private String getClassOfScore(int score) {

        if (0 < score && score < 10)
            return "0-10";
        else if (score < 20)
            return "10-20";
        else if (score < 30)
            return "20-20";
        else return "40-more";


    }

    void process(RowRepeatDB row, PrintWriter writer) throws StructureException {

        String[] units = row.getUnits().split(";");
        String pdbCode = row.getPdbCode();
        String pdbChain = row.getPdbChain();
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        Atom[] atoms = PdbTools.getAtomCAArray(structure
                .getChainByPDB(pdbChain));
        //System.out.println(">>++++" + pdbCode);

        Set<String> hashSet = this.tripeUnitsProcess(units, atoms, writer, row);

        for (String str : hashSet) {
            writer.write(str + "\n");
        }


//        MutilAlign mutilAlign = new MutilAlign();
//        for (int i = 0; i < units.length - 1; i++) {
//            try {
//                int starti = Integer.parseInt(units[i].split("-")[0]);
//                int endi = Integer.parseInt(units[i].split("-")[1]);
//
//                Atom[] atomi = this.getAtoms(atoms, starti, endi);
//
//				for (int j = i + 1; j < units.length; j++) {
//
//                    int startj = Integer.parseInt(units[j].split("-")[0]);
//                    int endj = Integer.parseInt(units[j].split("-")[1]);
//
//                    Atom[] atomj = this.getAtoms(atoms, startj, endj);
//
//                    AFPChain afpChain = mutilAlign.pairAlign(atomi, atomj);
//
//                    double rmsd = afpChain.getChainRmsd();
//                    int alnLeng = afpChain.getOptLength();
//                    int gapLeng = afpChain.getGapLen();
//                    int leng1 = afpChain.getCa1Length();
//                    int leng2 = afpChain.getCa2Length();
//                    int minL = Math.min(leng1, leng2);
//                    int maxL = Math.max(leng1, leng2);
//                    double aL = (double) alnLeng / maxL;
//                    double aS = (double) alnLeng / minL;
//
//                    double tmScore = afpChain.getTMScore();
//
//                    //write header
//
//                    // System.out.println(alnLeng + "\t" + df.format(rmsd) +
//                    // "\t"
//                    // + gapLeng + "\t" + leng1 + "\t" + leng2);
//
//                    writer.write(row.getEntry() + "(" + units[i] + " vs " + units[j] + ")" + "\t" + alnLeng + "\t" + df.format(rmsd) + "\t"
//                            + gapLeng + "\t" + leng1 + "\t" + leng2 + "\t" + minL + "\t" + maxL + "\t" + aL + "\t" + aS + "\t"
//                            + df.format(tmScore) + "\n");
//
//				}
//
//			} catch (Exception e) {
//				// TODO: handle exception
//			}
//
//		}

	}

	int getPosition(Atom[] atoms, int posNsq) {

		for (int i = 0; i < atoms.length; i++) {

			int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
			if (seqNumber == posNsq)
				return i;

		}

		return 0;
	}

	Atom[] getAtoms(Atom[] atoms, int startNsq, int endNsq) {

		List<Atom> subAtoms = new ArrayList<Atom>();

		for (Atom atom : atoms) {

			int seqNumber = atom.getGroup().getResidueNumber().getSeqNum();
			if (seqNumber >= startNsq && seqNumber <= endNsq)
				subAtoms.add(atom);

		}

		return subAtoms.toArray(new Atom[subAtoms.size()]);
	}
}
