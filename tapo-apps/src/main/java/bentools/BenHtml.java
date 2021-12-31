package bentools;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.cnrs.crbm.lib.analysis.ScoreAnalysis;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 7/2/2015.
 */
public class BenHtml {

    public static void main(String[] args) throws Exception {
        BenHtml benchmarking = new BenHtml();
        DataIO.writeToFile(benchmarking.generateHTML(), "F:\\App\\wamp\\www\\3DRepeat\\data\\report.html");
        //System.out.println(benchmarking.generateHTML());
    }


    public String generateHTML() throws Exception {

        StringBuilder builder = new StringBuilder();
        StringBuilder builderHeader = new StringBuilder();
        StringBuilder builderContent = new StringBuilder();
        Map<String, List<Repeat>> retrievedTRs = new ScoreAnalysis().readResult("data/benchmarkset/TRsRegion.o");
        //int a = 3;
        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        for (String row : releventTRs) {
            String[] data = row.split("\t");
            String pdb = data[0];

            //System.out.println(pdb);
            int startRegion = 0;
            int endRegion = 0;
            String region = data[2];
            String strClass = data[1];
            if (region.startsWith("-")) {
                region = region.substring(1, region.length());
                startRegion = (-1) * Integer.parseInt(region.split("-")[0]);
                endRegion = Integer.parseInt(region.split("-")[1]);

            } else {
                startRegion = Integer.parseInt(region.split("-")[0]);
                endRegion = Integer.parseInt(region.split("-")[1]);
            }
//            if (!pdb.equals("1ofl_A"))
//                continue;
            // pdb Code
            String pdbCode = pdb.substring(0, 4);
            String pdbChain = pdb.substring(5, 6);

            Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
            // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
            // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
            //
            Chain c = structure.getChainByPDB(pdbChain);
            //System.out.println(structure);
            Atom[] atoms = PdbTools.getAtomCAArray(c);

            startRegion = this.getPosition(atoms, startRegion);
            endRegion = this.getPosition(atoms, endRegion);
            if (endRegion < startRegion || endRegion == 0) {
                endRegion = atoms.length - 1;
            }


            int nRU = Integer.parseInt(data[4]);
            int L = (endRegion - startRegion + 1) / nRU;
            Repeat refRepeat = new Repeat();
            for (int i = 0; i < nRU; i++) {

                int start = i * L + startRegion;
                int end = i * L + startRegion + L - 1;
                if (end >= atoms.length)
                    end = atoms.length - 1;
                refRepeat.getRepeats().add(new RepeatContent(start, end));
            }


            /**
             * LOAD repeats which are predicted by TAPO
             */
            List<Repeat> tapoPredictedTRs = new ArrayList<Repeat>();
            if (retrievedTRs.containsKey(pdb)) {
                tapoPredictedTRs = retrievedTRs.get(pdb);
            }


            builderHeader.append("<div class=\"idDom\">" + pdb + " (Manually Annotated)</div>\n");
            for (Repeat repeat : tapoPredictedTRs) {
                // build one of repeat
                builderHeader.append("<div class=\"idDom\">" + pdb + " (predicted)</div>");

            }


            builderContent.append(buildOneRepeat(refRepeat, atoms, "ref") + "\n");
            builderContent.append("<div class=\"contDomClear\"></div>\n");
            for (Repeat repeat : tapoPredictedTRs) {
                // build one of repeat
                // check repeat

                //System.out.println(repeat.getFinderName());
                if (repeat.getFinderName().equals("CESymm"))// CE-Symm problem
                {
                    int nRU_Fix = repeat.getRepeats().size();
                    int L_Fix = (int) repeat.getAvgLength();
                    Repeat newRepeat = new Repeat();
                    for (int i = 0; i < nRU_Fix; i++) {
                        int start = i * L_Fix + repeat.getStart();
                        int end = i * L_Fix + repeat.getStart() + L_Fix - 1;
                        if (end < atoms.length - 1) {
                            newRepeat.getRepeats().add(new RepeatContent(start, end));
                        }

                    }
                    repeat = newRepeat;
                }


                builderContent.append(buildOneRepeat(repeat, atoms, "pre") + "\n");
                builderContent.append("<div class=\"contDomClear\"></div>\n");
            }


//            a--;
//            if (a < 0)
//                break;

        }
        builder.append("<div id=\"contIds\">");
        builder.append(builderHeader.toString());
        builder.append("</div>");
        builder.append("<div  id=\"contSeqDomGlobal\">");
        builder.append(builderContent.toString());
        builder.append("</div>");
        return builder.toString();

    }

    private String buildOneRepeat(Repeat repeat, Atom[] atoms, String type) {
        StringBuffer buffer = new StringBuffer();
        buffer.append("<div class=\"contDom\">\n");
        List<RepeatContent> listTR = repeat.getRepeats();
        buffer.append(this.buildLine(convertToNewSize(listTR.get(0).getStart(), atoms), atoms));
        double size = this.convertToNewSize(listTR.get(0).getStart(), atoms);
        double accumunatedPoint = size;
        int index = 0;
        for (int i = 0; i < listTR.size() - 1; i++) {
            //build block
            // space between 2 block
            size = this.convertToNewSize(listTR.get(i).size(), atoms);
            int space = listTR.get(i + 1).getStart() - listTR.get(i).getEnd() - 1;
            accumunatedPoint += (size + convertToNewSize(space, atoms));
            if (i % 2 == 1)
                type = "pre";
            else
                type = "ref";
            buffer.append(buildRepeatBlock(listTR.get(i), atoms, type));
            buffer.append(buildLine(convertToNewSize(space, atoms), atoms));
            index++;

        }
        if (index % 2 == 1)
            type = "pre";
        else
            type = "ref";
        size = this.convertToNewSize(listTR.get(listTR.size() - 1).size(), atoms);
        accumunatedPoint += size;
        buffer.append(buildRepeatBlock(listTR.get(listTR.size() - 1), atoms, type));
        buffer.append(this.buildLine(SCALE_FACTOR - accumunatedPoint, atoms));
        //buffer.append(this.buildLine(convertToNewSize(atoms.length - listTR.get(listTR.size() - 1).getEnd(), atoms), atoms));
//        if (accumunatedPoint < SCALE_FACTOR -size) {
//            buffer.append(buildRepeatBlock(listTR.get(listTR.size() - 1), atoms, type));
//            buffer.append(this.buildLine(SCALE_FACTOR - accumunatedPoint, atoms));
//        } else {
//            buffer.append(this.buildLine(SCALE_FACTOR - accumunatedPoint, atoms));
//        }

        buffer.append("</div>\n");
        return buffer.toString();
    }

    private String buildRepeatBlock(RepeatContent r, Atom[] atoms, String type) {
        try {
            double size = this.convertToNewSize(r.size(), atoms);
            int start = getResSeq(r.getStart(), atoms);
            int end = getResSeq(r.getEnd(), atoms);

            String css = "seqDom";
            if (type.equals("ref"))
                css = "seqDomRef";

            return "<div  class=\"" + css + "\" title=\"" + start + " - " + end + "\"  style=\"width: " + size + "px;\"></div>\n";
        } catch (Exception ex) {
            return "";
        }
    }

    private String buildLine(double size, Atom[] atoms) {

        return "<div  class=\"basicSeqDom\" style=\"width: " + size + "px;\"></div>\n";
    }

    private int convertToNewSize(int oldSize, Atom[] atoms) {
        double scale = (double) SCALE_FACTOR / atoms.length;
        return (int) Math.floor(oldSize * scale);
//        //if (newSize <= 6) newSize = 1;
//        return newSize;
    }

    private final static int SCALE_FACTOR = 600;

    private static int getResSeq(int pos, Atom[] atoms) {

        return atoms[pos].getGroup().getResidueNumber().getSeqNum();
    }

    int getPosition(Atom[] atoms, int posNsq) {
        for (int i = 0; i < atoms.length; i++) {

            int seqNumber = atoms[i].getGroup().getResidueNumber().getSeqNum();
            if (seqNumber == posNsq)
                return i;

        }
        return 0;
    }

}
