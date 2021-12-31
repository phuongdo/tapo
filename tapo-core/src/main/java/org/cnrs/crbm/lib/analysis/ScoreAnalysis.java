package org.cnrs.crbm.lib.analysis;

import org.cnrs.crbm.lib.conf.ThresholdConfig;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ProteinCSVReader;
import org.cnrs.crbm.lib.io.Row;
import org.cnrs.crbm.lib.repeats.clusters.ClusterLocation;
import org.cnrs.crbm.lib.repeats.clusters.ClusterRepeat;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.trsfinder.RepeatContent;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 2/20/2015.
 */
public class ScoreAnalysis {

    public static void main(String[] args) throws Exception {


        //new ScoreAnalysis().run();


    }


    public Repeat convertToRepeat(String line) {
        Repeat repeat = new Repeat();
        String[] data = line.split("\t");
        try {

            String[] units = data[6].split(";");
            repeat.setScore(Double.parseDouble(data[7]));
            repeat.setRankScore(Double.parseDouble(data[8]));
            repeat.setFinderName(data[1]);


            for (String unit : units) {
                int start = 0;
                int end = 0;
                if (unit.startsWith("-")) {
                    unit = unit.substring(1, unit.length());
                    start = (-1) * Integer.parseInt(unit.split("-")[0]);
                    end = Integer.parseInt(unit.split("-")[1]);
                    repeat.getRepeats().add(new RepeatContent(start, end));

                } else {
                    start = Integer.parseInt(unit.split("-")[0]);
                    end = Integer.parseInt(unit.split("-")[1]);
                    repeat.getRepeats().add(new RepeatContent(start, end));
                }
            }

        } catch (Exception ex) {

        }

        return repeat;

    }

    //static double QA_THRES = 0.6;

    public Map<String, List<Repeat>> readResult(String fileDir) throws Exception {

        //String fileDir = "data/tapo/ExperimentalResult.o";
        List<String> lines = DataIO.readLines(fileDir);
        Map<String, List<Repeat>> map = new HashMap<String, List<Repeat>>();
        for (String line : lines) {
            //System.out.println(line);
            if (line.length() < 10|| line.contains("No-TRs"))
                continue;
            try {
                String[] data = line.split("\t");
                String pdbEntry = data[0];
//                if (pdbEntry.equals("3ams_A"))
//                    System.out.println();
                Repeat repeat = this.convertToRepeat(line);
                if (repeat.getScore() >= 0.0) {
                    if (map.containsKey(pdbEntry)) {
                        List<Repeat> l = map.get(pdbEntry);
                        l.add(repeat);
                        map.put(pdbEntry, l);
                    } else {
                        List<Repeat> l = new ArrayList<Repeat>();
                        l.add(repeat);
                        map.put(pdbEntry, l);
                    }
                }

            } catch (Exception ex) {
                System.out.println();
            }
        }
        return map;

    }

    public void run() throws Exception {

        Map<String, List<Repeat>> map = this.readResult("data/tapo/ExperimentalResult.o");
        // read from francois
        ProteinCSVReader csvReader = new ProteinCSVReader();
        List<Row> rowsFrancois = csvReader.getData("data/tapo/FrancoisData1802.csv");
        for (Row row : rowsFrancois) {
            if (map.containsKey(row.getProtein().trim())) {


//                if(row.getProtein().equals("1h3i_A"))
//                    System.out.println();
                List<Repeat> list = map.get(row.getProtein().trim());
                List<Repeat> allTRs = new ArrayList<Repeat>();

                for (Repeat repeat : list) {
                    if (repeat.getScore() >= ThresholdConfig.QA_THRES)
                        allTRs.add(repeat);
                }
                if (allTRs.size() > 0) {
                    List<ClusterRepeat> clusters = new ClusterLocation().cluster(allTRs);
                    ClusterLocation clusterLocation = new ClusterLocation();

                    for (ClusterRepeat c : clusters) {
                        StringBuffer contentBuilder = new StringBuffer();
                        // System.out.println("-----");
                        c.sortScoreDESC();
                        //get top
                        Repeat top = c.getRepeats().get(0);
                        int start = top.getStart();
                        int end = top.getEnd();
                        // pdbId_chain Finder avgLeng nUnits RL
                        StringBuffer unitsBuffer = new StringBuffer();
                        for (RepeatContent unit : top.getRepeats()) {
                            unitsBuffer.append(unit.getStart() + "-" + unit.getEnd() + ";");
                        }
                        String unitsStr = unitsBuffer.toString();
                        unitsStr = unitsStr.substring(0, unitsStr.length() - 1);
                        contentBuilder.append(
                                +top.getRepeats().size() + "\t"
                                        + NumberFormatUtils.format(top.getAvgLength()) + "\t" + start
                                        + "\t" + end + "\t" + unitsStr + "\t" + NumberFormatUtils.format(top.getScore()) + "\t" + NumberFormatUtils.format(top.getRankScore()))
                        ;


//                        Repeat repeatXX = new Repeat();
//                        repeatXX.getRepeats().add(new RepeatContent(row.getPdbBegin(), row.getPdbEnd()));
                        double overlap = clusterLocation.overlap(new Region(row.getPdbBegin(),row.getPdbEnd()), new Region(top.getStart(),top.getEnd()));
                        System.out.println(row.getProtein() + "\t" + row.getByEyeTot() + "\t" + row.getPdbBegin() + "\t" + row.getPdbEnd() + "\t" + contentBuilder.toString() + "\t" + NumberFormatUtils.format(overlap));


                    }


                }
            } else {
                System.out.println(row.getProtein() + "\t" + row.getByEyeTot() + "\t" + row.getPdbBegin() + "\t" + row.getPdbEnd() + "\t0\t0\t0\t0\t0\t0\t0\t0");

            }

        }

    }


}
