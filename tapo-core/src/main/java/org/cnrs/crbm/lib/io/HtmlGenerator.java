package org.cnrs.crbm.lib.io;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.cnrs.crbm.lib.conf.Dir;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.ColorUtils;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 4/17/2015.
 */
public class HtmlGenerator {

    String pdbCode = "";
    String pdbChain = "";
    String pid = "";
    Atom[] atoms = null;
    private final static int SCALE_FACTOR = 600;

    public void setSHOW_OPT(String SHOW_OPT) {
        this.SHOW_OPT = SHOW_OPT;
    }

    private String SHOW_OPT = "All";


    StringBuffer buffer = new StringBuffer();


    public static void main(String[] args) {
        HtmlGenerator generator = new HtmlGenerator("1c1g", "C", "55352373a5cc3");

        generator.buildHTML();
        //DataIO.writeToFile("output/demo.html", generator.html());
        System.out.println(generator.html());

    }


    public HtmlGenerator(String pdbCode, String pdbChain, String pid) {
        this.pdbCode = pdbCode;
        this.pdbChain = pdbChain;
        this.pid = pid;
        // version 1.0.2 . support get structure directly from local install of
        // PDB database.
        Structure structure = PdbTools.getStructureFromLocalPdb(pdbCode);
        // Atom[] cbAtoms = PdbTools.getAtomCBArray(StructureTools
        // .getAtomCAArray(structure.getChainByPDB(pdbChain)));
        //

        try {

            Chain c = structure.getChainByPDB(pdbChain);
            //System.out.println(structure);
            atoms = PdbTools.getAtomCAArray(c);
        } catch (StructureException e) {
            e.printStackTrace();
        }


    }

    ReadFasta readFasta = new ReadFasta();

    private List<Repeat> filter(List<Repeat> repeats) {

        List<Repeat> repeatsFilter = new ArrayList<Repeat>();

        if (this.SHOW_OPT.equals("All"))
            return repeats;
        else {

            for (Repeat repeat : repeats) {
                if (repeat.getCluster().contains("selected")) {
                    repeatsFilter.add(repeat);
                }
            }
            return repeatsFilter;
        }

    }

    public void buildHTML() {

        String dirOutput = Dir.TMP_DIR + "/" + pid + "_" + pdbCode + "" + pdbChain + ".o";
        OutputReader reader = new OutputReader();

        try {
            Map<String, TaPoFastaFormat> fastaTR = readFasta.readTaPoFastaFormat(dirOutput);
            //Map<String, List<Repeat>> pdbOut = reader.readResult(dirOutput);
            TaPoFastaFormat taPoFastaFormat = fastaTR.get(pdbCode + "_" + pdbChain);
            //List<Repeat> repeats = pdbOut.get(pdbCode + "_" + pdbChain);
            List<Repeat> repeats = taPoFastaFormat.getRepeats();
            repeats = this.filter(repeats);
            //

            // build description.
            if (!this.SHOW_OPT.equals("All")) {
                String status = "GLOBULAR";
                String statusColor = "red";
                if(taPoFastaFormat.is3DRepeat()) {
                    status = "REPEAT";
                    statusColor = "green";
                }
                buffer.append("<div>");
                buffer.append("Query structure: <b>" + pdbCode + "</b> Chain : <b>" + pdbChain + "</b> <br /> " +
                        "SVM-Score: <b>" + taPoFastaFormat.getSvmScore() + "</b><br/>\n" +
                        "Others: " + taPoFastaFormat.getCombineScore() + "<br/>\n" +
                        "Status:<font\n" +
                        "style=\"color: "+statusColor+"\"> <b>"+status+"</b></font><br/><br/>\n");

                buffer.append("</div>");
            }

            buffer.append("<div id=\"contIds\">");
            for (Repeat repeat : repeats) {
                // build one of repeat
                buffer.append("<div class=\"idDom\">" + buildViewMSA(repeat) + "</div>");

            }
            buffer.append("</div>");
            buffer.append("<div  id=\"contSeqDomGlobal\">");
            for (Repeat repeat : repeats) {
                // build one of repeat
                buffer.append(buildOneRepeat(repeat) + "\n");
                buffer.append("<div class=\"contDomClear\"></div>\n");
            }
            buffer.append("</div>");

        } catch (Exception e) {
//            e.printStackTrace();
        }


    }

    private String buildOneRepeat(Repeat repeat) {
        StringBuffer buffer = new StringBuffer();
        buffer.append("<div class=\"contDom\">\n");
        List<RepeatContent> listTR = repeat.getRepeats();
        int size = this.convertToNewSize(listTR.get(0).getStart());
        buffer.append(this.buildLine(size));
        double accumunatedPoint = size;
        int index = 0;
        String type = "pre";
        for (int i = 0; i < listTR.size() - 1; i++) {
            size = this.convertToNewSize(listTR.get(i).size());
            int space = listTR.get(i + 1).getStart() - listTR.get(i).getEnd() - 1;
            accumunatedPoint += (size + convertToNewSize(space));
            //build block
            if (i % 2 == 0)
                type = "pre";
            else
                type = "ref";
            buffer.append(buildRepeatBlock(listTR.get(i), size, type));
            // space between 2 block
            buffer.append(buildLine(convertToNewSize(space)));
            index++;

        }
        if (index % 2 == 0)
            type = "pre";
        else
            type = "ref";
        size = this.convertToNewSize(listTR.get(listTR.size() - 1).size());
        accumunatedPoint += size;
        buffer.append(buildRepeatBlock(listTR.get(listTR.size() - 1), size, type));
        buffer.append(this.buildLine(SCALE_FACTOR - accumunatedPoint));
        buffer.append("</div>\n");
        return buffer.toString();
    }
//    private String buildOneRepeat(Repeat repeat) {
//        StringBuffer buffer = new StringBuffer();
//        buffer.append("<div class=\"contDom\">\n");
//
//        List<RepeatContent> listTR = repeat.getRepeats();
//
//        buffer.append(this.buildLine(convertToNewSize(listTR.get(0).getStart())));
//
//        for (int i = 0; i < listTR.size() - 1; i++) {
//
//            //build block
//            buffer.append(buildRepeatBlock(listTR.get(i)));
//
//            // space between 2 block
//            int space = listTR.get(i + 1).getStart() - listTR.get(i).getEnd();
//            buffer.append(buildLine(convertToNewSize(space)));
//        }
//        buffer.append(buildRepeatBlock(listTR.get(listTR.size() - 1)));
//        buffer.append(this.buildLine(convertToNewSize(atoms.length - listTR.get(listTR.size() - 1).getEnd())));
//        buffer.append("</div>\n");
//        return buffer.toString();
//    }

    private String buildLine(double size) {

        return "<div  class=\"basicSeqDom\" style=\"width: " + size + "px;\"></div>\n";
    }

    private String buildRepeatBlock(RepeatContent r, int size, String type) {

        //double size = this.convertToNewSize(r.size());
        int start = getResSeq(r.getStart(), atoms);
        int end = getResSeq(r.getEnd(), atoms);
        String css = "seqDom";
        if (type.equals("ref"))
            css = "seqDomRef";
        return "<div  class=\"" + css + "\" title=\"" + start + " - " + end + "\"  style=\"width: " + size + "px;\"></div>\n";
        // return "<div  class=\"seqDom\" title=\"" + start + " - " + end + "\"  style=\"width: " + size + "px;\"></div>\n";
    }

    public String buildViewMSA(Repeat repeat) {
        String jmolTemp = "select *;color gray;";
        StringBuilder builderJmol = new StringBuilder();
        StringBuilder builderMutilAlign = new StringBuilder();
        StringBuilder builderContent = new StringBuilder();

        int old = 0;
        for (RepeatContent r : repeat.getRepeats()) {

            int start = getResSeq(r.getStart(), atoms);
            int end = getResSeq(r.getEnd(), atoms);
            builderJmol.append("select " + start + "-" + end + ";");
            builderMutilAlign.append(r.getStart() + "-"
                    + r.getEnd() + ";");

            int size = ColorUtils.getInstance().getColors().size();

            int abudant = old % size;

            int red = ColorUtils.getInstance().getColors()
                    .get(abudant).getRed();
            int green = ColorUtils.getInstance().getColors()
                    .get(abudant).getGreen();
            int blue = ColorUtils.getInstance().getColors()
                    .get(abudant).getBlue();

            builderJmol.append("color [" + red + "," + green + ","
                    + blue + "];");


            old++;
        }

//      class='viewIcon'
//        String linkView = "<a href=\"javascript:callJSmolScript('"
//                + jmolTemp
//                + builderJmol.toString()
//                + "');view_mutilalign_ajax('"
//                + builderMutilAlign.toString().substring(0,
//                builderMutilAlign.length() - 1)
//                + "');\" >View</a>\n";

        String linkView = "<a href=\"javascript:visualizeRepeats('"
                + jmolTemp
                + builderJmol.toString()
                + "','"
                + builderMutilAlign.toString().substring(0,
                builderMutilAlign.length() - 1)
                + "');\" >" + repeat.getFinderName() + "</a>\n";

        // summary infor

        String infor = " [Q:" + NumberFormatUtils.format(repeat.getScore()) + " R:" + NumberFormatUtils.format(repeat.getRankScore()) + "]";
        //String infor = " [" + repeat.getFinderName() + " Q-score:" + repeat.getScore() + "]";
        return infor + " " + linkView;
    }


    private int convertToNewSize(int oldSize) {
        double scale = (double) SCALE_FACTOR / atoms.length;
        return (int) Math.round(oldSize * scale);
    }

    public String html() {
        return buffer.toString();
    }

    private static int getResSeq(int pos, Atom[] atoms) {

        return atoms[pos].getGroup().getResidueNumber().getSeqNum();
    }


}
