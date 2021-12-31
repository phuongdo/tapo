package org.cnrs.crbm.lib.utils;

import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.repeats.clusters.ClusterRepeat;
import org.cnrs.crbm.lib.repeats.shortTRs.LocationMatch;
import org.cnrs.crbm.lib.trsfinder.Finder;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;

import java.util.List;

public class HtmlUtils {


    /**
     * html generated for all tandem repeats
     */
    public static String getOutHtmlForTRs(List<ClusterRepeat> clusters, List<LocationMatch> listFromDetector, Atom[] atoms) {

        String jmolTemp = "select *;color gray;";
        StringBuilder builder = new StringBuilder();

        builder.append("<tr><td>RL</td><td>cluster</td><td>NoTRs</td><td>avgLength</td><td>Location</td><td>RUs</td><td>QA score</td><td>R score</td><td></td></tr>");
        try {
            int oldX = 0;
            //for short detection.

            for (LocationMatch match : listFromDetector) {
                StringBuilder builderJmolX = new StringBuilder();
                //StringBuilder builderMutilAlignX = new StringBuilder();
                StringBuilder builderContentX = new StringBuilder();
                int start = match.getStart();
                int end = match.getEnd();
                int seqStart = atoms[start].getGroup()
                        .getResidueNumber().getSeqNum();
                int seqEnd = atoms[end].getGroup()
                        .getResidueNumber().getSeqNum();

                builder.append("<tr>");

                int startR = getResSeq(start, atoms);
                int endR = getResSeq(end, atoms);
                builderJmolX.append("select " + startR + "-" + endR + ";");
//                builderMutilAlignX.append(start + "-"
//                        + end + ";");

                int size = ColorUtils.getInstance().getColors().size();

                int abudant = oldX % size;

                int red = ColorUtils.getInstance().getColors()
                        .get(abudant).getRed();
                int green = ColorUtils.getInstance().getColors()
                        .get(abudant).getGreen();
                int blue = ColorUtils.getInstance().getColors()
                        .get(abudant).getBlue();

                builderJmolX.append("color [" + red + "," + green + ","
                        + blue + "];");

                builderContentX.append("<td>" + match.getMeg() + "</td><td>clusX</td>" + "<td>"
                        + 1 + "</td><td>"
                        + (NumberFormatUtils.format(match.getEnd() - match.getStart() + 1)) + "</td><td>" + seqStart
                        + "-" + seqEnd + "</td><td>" + seqStart + "-" + seqEnd + "</td><td>" + NumberFormatUtils.format(0.0) + "</td><td>" + NumberFormatUtils.format(0.0) + "</td>");
                String linkView = "<a href=\"javascript:document.jmol.script('"
                        + jmolTemp
                        + builderJmolX.toString()
                        + "');\" ><font color='red'>  View "

                        + "</font></a>";

                builder.append(builderContentX.toString() + "<td>" + linkView + "</td>");
                ;

                oldX++;
                builder.append("</tr>");
            }

            int indexG = 1;
            for (ClusterRepeat c : clusters) {
                // System.out.println("-----");
                c.sortScoreDESC();
                //get top
                Repeat top = c.getRepeats().get(0);
                Repeat repeat = top;
                int index = 0;
                //for (Repeat repeat : c.getRepeats()) {
                builder.append("<tr>");


                int startR = repeat.getStart();
                int endR = repeat.getEnd();
                // pdbId_chain Finder avgLeng nUnits RL
                int seqStart = atoms[startR].getGroup()
                        .getResidueNumber().getSeqNum();
                int seqEnd = atoms[endR].getGroup()
                        .getResidueNumber().getSeqNum();
                StringBuffer unitsBuffer = new StringBuffer();
                for (RepeatContent unit : repeat.getRepeats()) {

                    int seqS = atoms[unit.getStart()].getGroup()
                            .getResidueNumber().getSeqNum();
                    int seqE = atoms[unit.getEnd()].getGroup()
                            .getResidueNumber().getSeqNum();

                    unitsBuffer.append(seqS + "-" + seqE + ";");
                }
                String unitsStr = unitsBuffer.toString();
                unitsStr = unitsStr.substring(0, unitsStr.length() - 1);

                try {
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

                    builderContent.append("<td>" + "tapo" + "</td><td>" + "clus" + indexG
                            + "</td><td>"
                            + repeat.getRepeats().size() + "</td><td>"
                            + NumberFormatUtils.format(repeat.getAvgLength()) + "</td><td>" + seqStart
                            + "-" + seqEnd + "</td><td>" + unitsStr + "</td><td>" + NumberFormatUtils.format(repeat.getScore()) + "</td><td>" + NumberFormatUtils.format(repeat.getRankScore()) + "</td>")
                    ;


                    String linkView = "<a href=\"javascript:document.jmol.script('"
                            + jmolTemp
                            + builderJmol.toString()
                            + "');view_mutilalign_ajax('"
                            + builderMutilAlign.toString().substring(0,
                            builderMutilAlign.length() - 1)
                            + "');\" ><font color='red'>  View "

                            + "</font></a>\n";

                    builder.append(builderContent.toString() + "<td>" + linkView + "</td>");
                    builder.append("</tr>");
                    index++;

                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }


                //}

                indexG++;
            }

        } catch (Exception ex) {

        }


        return builder.toString();


    }


    /**
     * html generated for each finder
     *
     * @param finder
     * @param program
     * @return
     */
    public static String getOutHtml(Finder finder, String program, Atom[] atoms) {
        String jmolTemp = "select *;color gray;";

        List<Repeat> repeats = finder.getRepeats();
        StringBuilder builder = new StringBuilder();
        builder.append("<td valign='top'>");
        builder.append("<b>" + program.toUpperCase() + "</b><br/>\n");
        int type = 1;
        if (repeats.size() > 0 && finder.isTRs()) {

            for (Repeat repeat : repeats) {

                try {
                    StringBuilder builderJmol = new StringBuilder();
                    StringBuilder builderMutilAlign = new StringBuilder();
                    StringBuilder builderContent = new StringBuilder();
                    builderContent.append("avgL:"
                            + NumberFormatUtils.format(repeat.getAvgLength())
                            + "<br/>\n");
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

                        // if (abudant == 0)
                        // builderJmol.append("color yellow;");
                        // else if (abudant == 1)
                        // builderJmol.append("color green;");
                        // else if (abudant == 2)
                        // builderJmol.append("color red;");
                        // else if (abudant == 3)
                        // builderJmol.append("color blue;");
                        // else if (abudant == 4)
                        // builderJmol.append("color orange;");

                        builderContent.append(""
                                + getResSeq(r.getStart(), atoms) + "-"
                                + getResSeq(r.getEnd(), atoms) + "<br/>\n");

                        old++;
                    }

                    String linkView = "<a href=\"javascript:document.jmol.script('"
                            + jmolTemp
                            + builderJmol.toString()
                            + "');view_mutilalign_ajax('"
                            + builderMutilAlign.toString().substring(0,
                            builderMutilAlign.length() - 1)
                            + "');\" ><font color='red'> >RL: "
                            + type
                            + "</font></a><br/>\n";

                    builder.append(linkView + builderContent.toString());
                    type++;
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
        } else {
            builder.append("NO<br/>\n");
        }
        builder.append("</td>");
        builder.append("<td></<td>");
        return builder.toString();

    }

    /**
     * html genderate for combined methods
     *
     * @param finder
     * @param program
     * @return
     */
    public static String getOutHtmlForCombinedMethod(Finder finder,
                                                     String program, Atom[] atoms) {

        int nrAtoms = atoms.length;

        List<Repeat> repeats = finder.getRepeats();
        StringBuilder builder = new StringBuilder();
        if (repeats.size() > 0 && finder.isTRs()) {

            int index = 1;
            for (Repeat repeat : repeats) {
                try {
                    // for each of repeat type
                    // convert region?

                    double first = ((double) repeat.getRepeats().get(0)
                            .getStart() / nrAtoms) * 100;
                    double second = ((double) repeat.getAvgLength()
                            * repeat.getRepeats().size() / nrAtoms) * 100;
                    double third = 100 - first - second;

                    builder.append("<tr>\n");
                    builder.append("<td>" + program + "</td>\n");
                    builder.append("<td>Type " + index + "</td>\n");
                    builder.append("<td class='graph'><table  width='100%'>\n");
                    builder.append("<tbody>\n");
                    builder.append("<tr>\n");
                    builder.append("<td class='add' style='width: " + first
                            + "%;'></td>\n");
                    builder.append("<td class='rem' style='width: " + second
                            + "%;'></td>\n");
                    builder.append("<td class='add' style='width: " + third
                            + "%;'></td>\n");
                    builder.append("</tr>\n");
                    builder.append("</tbody>\n");
                    builder.append("</table></td>\n");
                    builder.append("</tr>\n");

                    index++;

                } catch (Exception e) {
                    System.out.println(e.getMessage());
                }
            }
        } else {
            // no repeat is found here

        }
        return builder.toString();

    }

    private static int getResSeq(int pos, Atom[] atoms) {

        return atoms[pos].getGroup().getResidueNumber().getSeqNum();
    }

}
