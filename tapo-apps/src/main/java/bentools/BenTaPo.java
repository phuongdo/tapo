package bentools;

import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.repeats.CombineScore;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.repeats.clusters.ClusterLocation;
import org.cnrs.crbm.lib.trsfinder.Region;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 9/11/2015.
 */
public class BenTaPo {
    public static void main(String[] args) throws Exception {


        new BenTaPo().benchmarkRepeatLength();
//        new BenTaPo().benchmark();


//        for (Map.Entry<String, TaPoFastaFormat> entry : fastaTRs.entrySet()) {
//            TaPoFastaFormat taPoFastaFormat = entry.getValue();
//
//            taPoFastaFormat.toString().
//        }

    }


    public void benchmarkRepeatLength() throws Exception {

        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/benchmarkset/tapov1.1.1.output");
        ClusterLocation clusterLocation = new ClusterLocation();
        List<String> releventTRs = DataIO.readLines("data/benchmarkset/TRsRegion.ab.in");
        for (String row : releventTRs) {
            try {
                String[] data = row.split("\t");
                String pdb = data[0];
                String pdbCode = pdb.substring(0, 4);
                String pdbChain = pdb.substring(5, 6);

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

                Region refRegion = new Region(startRegion, endRegion);
                double avgRefL = Double.parseDouble(data[3]);


                StringBuffer buffer = new StringBuffer();

                String strStart = pdb + "\t" + refRegion + "\t" + avgRefL + "\t";
                String falsePositive = "0.00\t0-0\t0.00";
                if (fastaTRs.containsKey(pdb)) {

                    TaPoFastaFormat taPoFastaFormat = fastaTRs.get(pdb);

                    if (taPoFastaFormat.is3DRepeat()) {
                        List<Repeat> lstRepeat = new ArrayList<Repeat>();
                        for (Repeat repeat : taPoFastaFormat.getRepeats()) {
                            if (repeat.getCluster().contains("selected")) {
                                lstRepeat.add(repeat);
                            }
                        }

                        // get

                        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                        Atom[] atoms = repeatFinder.getAtoms();

                        String display = "";
                        boolean found = false;
                        for (Repeat repeat : lstRepeat) {

                            int start = PdbTools.getResSeq(atoms, repeat.getStart());
                            int end = PdbTools.getResSeq(atoms, repeat.getEnd());
                            Region predRegion = new Region(start, end);
                            double cover = clusterLocation.cover(predRegion, refRegion);
                            if (cover >= 0.5) {
                                // true positive
                                //buffer.append(NumberFormatUtils.format(cover) + "\t" + start + "-" + end + "\t" + NumberFormatUtils.format(repeat.getAvgLength()));
                                System.out.println(strStart + "1\t" + NumberFormatUtils.format(cover) + "\t" + start + "-" + end + "\t" + NumberFormatUtils.format(repeat.getAvgLength()));
                                found = true;

                            }

                        }

                        if (!found) {
                            for (Repeat repeat : lstRepeat) {
                                int start = PdbTools.getResSeq(atoms, repeat.getStart());
                                int end = PdbTools.getResSeq(atoms, repeat.getEnd());
                                Region predRegion = new Region(start, end);
                                double cover = clusterLocation.cover(predRegion, refRegion);
                                System.out.println(strStart + "0\t" + NumberFormatUtils.format(cover) + "\t" + start + "-" + end + "\t" + NumberFormatUtils.format(repeat.getAvgLength()));
                            }

                        }

                    } else {
//                        System.out.println(strStart + "0\t" + falsePositive);
                    }


                } else {
//                    System.out.println(strStart + falsePositive);
                }


//                System.out.println(pdb + "\t" + refRegion + "\t" + avgRefL + "\t" + buffer.toString());


            } catch (Exception ex) {
//                ex.printStackTrace();
            }
        }

    }

    public void benchmark() throws Exception {

        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/benchmarkset/tapov1.1.1.output");
        List<String> datasetTRs = DataIO.readLines("data/benchmarkset/pub/Dataset-TR.txt");
        // print header
        System.out.println("protein,3DTRs," + CombineScore.toHeader());
        for (String row : datasetTRs) {
            try {
                TaPoFastaFormat taPoFastaFormat = fastaTRs.get(row);
                System.out.println(row + ",1," + taPoFastaFormat.getCombineScore().toCSVFormat());
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        List<String> datasetNonTRs = DataIO.readLines("data/benchmarkset/pub/Dataset-non-TR-v.1.1.0.txt");
        for (String row : datasetNonTRs) {
            try {
                TaPoFastaFormat taPoFastaFormat = fastaTRs.get(row);
                System.out.println(row + ",0," + taPoFastaFormat.getCombineScore().toCSVFormat());
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }

    }


}
