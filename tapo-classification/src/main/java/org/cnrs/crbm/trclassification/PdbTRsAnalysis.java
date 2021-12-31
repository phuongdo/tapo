package org.cnrs.crbm.trclassification;


import org.biojava.nbio.structure.Atom;
import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.repeats.RepeatFinder;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.utils.NumberFormatUtils;
import org.cnrs.crbm.lib.utils.PdbTools;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.util.*;

/**
 * Created by pdoviet on 7/20/2015.
 */
public class PdbTRsAnalysis {

    public static void main(String[] args) throws Exception {

        new PdbTRsAnalysis().classifyTRsBasedOnCATH();
//        new PdbTRsAnalysis().classifyTRsBasedOnSCOPe();
//        new PdbTRsAnalysis().processListOfStructure();
        //System.out.println(new PdbTRsAnalysis().parseDomain("a.1.1.1"));


    }

    public String convertToPymolVisualizationFormat(Repeat repeat, String pdbId) {
        StringBuffer buffer = new StringBuffer();
        String pdbCode = pdbId.substring(0, 4);
        String pdbChain = pdbId.substring(5, 6);
        //extract atoms
        RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
        repeatFinder.setOutputRelativeScore(false);
        return repeatFinder.buildOutputOneRepeat(repeat).toString();
    }

    private void processListOfStructure() throws Exception {
        Set<String> processList = this.getProcessList();
        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/resultLargeScale/tapo-v1.1.2/tapo_PDBTRs_June2015_v1.1.2.out");
        StringBuilder builder = new StringBuilder();
        //StringBuilder builder = new StringBuilder();

        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = fastaTRs.size();
        for (Map.Entry<String, TaPoFastaFormat> entry : fastaTRs.entrySet()) {
//            System.out.println(entry.getKey());
            //System.out.println(">>>>>>" + entry.getValue());
            bar.update(process, sizeOfProcess);
            process++;

            String pdbId = entry.getKey();
            String pdbCode = pdbId.substring(0, 4);
            String pdbChain = pdbId.substring(5, 6);
            if (!processList.contains(pdbId))
                continue;
            try {
                TaPoFastaFormat taPoFastaFormat = entry.getValue();
                if (taPoFastaFormat.is3DRepeat()) {
                    if (taPoFastaFormat.getRepeats().size() > 0) {
                        Repeat repeat = taPoFastaFormat.getRepeats().get(0);
//                        builder.append(this.convertToPymolVisualizationFormat(repeat, pdbId) + "\n");
                        ProteinSCOP scopeDomain = ScopeClassifierML.getInstance().classifyPdb(pdbCode, pdbChain, repeat.getStart(), repeat.getEnd());
                        builder.append(pdbId + "\t" + repeat.getStart() + "-" + repeat.getEnd() + "\t" + scopeDomain + "\t" + repeat.getAvgLength() + "\n")
                        ;
                    }
                }
            } catch (Exception ex) {
            }

        }//end


        DataIO.writeToFile(builder.toString(), "output/TrsUnkPDB.SCOP.txt");


    }

    public Set<String> getProcessList() {

        String fileDir = "data/trclassification/UnclassifyProtein";
        List<String> lines = DataIO.readLines(fileDir);
        Set<String> set = new HashSet<String>();
        for (String line : lines) {
            set.add(line);
        }
        return set;
    }

    public void classifyTRsBasedOnSCOPe() throws Exception {

//        PredLenMod predLenMod = new PredLenMod();
        // read tapo large-scale output
        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/resultLargeScale/tapo-v1.1.2/tapo_PDBTRs_June2015_v1.1.2.out");
        StringBuilder builder = new StringBuilder();
        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = fastaTRs.size();
        for (Map.Entry<String, TaPoFastaFormat> entry : fastaTRs.entrySet()) {
//            System.out.println(entry.getKey());
            //System.out.println(">>>>>>" + entry.getValue());
            bar.update(process, sizeOfProcess);
            process++;
            try {
                String pdbId = entry.getKey();
                String pdbCode = pdbId.substring(0, 4);
                String pdbChain = pdbId.substring(5, 6);
                TaPoFastaFormat taPoFastaFormat = entry.getValue();
                if (taPoFastaFormat.is3DRepeat()) {
                    if (taPoFastaFormat.getRepeats().size() > 0) {

//                        List<RepeatRegion> regions = predLenMod.predictByTAPOFormat(taPoFastaFormat);
//                        for (RepeatRegion region : regions) {
//                            //Repeat repeat = taPoFastaFormat.getRepeats().get(0);
//                            String scopeDomain = ScopeClassifierML.getInstance().classifyPdb(pdbCode, pdbChain, region.getStart(), region.getEnd());
////                            System.out.println(pdbId + "\t" + region.getSeqStart() + "-" + region.getSeqEnd() + "\t" + parseDomain(scopeDomain) + "\t" + region.getPredLen() + "\t" + region.getPredRU());
//
//                            builder.append(pdbId + "\t" + region.getSeqStart() + "-" + region.getSeqEnd() + "\t" + parseDomain(scopeDomain) + "\t" + region.getPredLen() + "\t" + region.getPredRU() + "\n");
//                        }

                        //  analysis the seleteted tandem repeat candidates
                        List<Repeat> lstRepeat = new ArrayList<Repeat>();
                        for (Repeat repeat : taPoFastaFormat.getRepeats()) {
                            if (repeat.getCluster().contains("selected")) {
                                lstRepeat.add(repeat);
                            }
                        }
                        for (Repeat repeat : lstRepeat) {
                            ProteinSCOP scopeDomain = ScopeClassifierML.getInstance().classifyPdb(pdbCode, pdbChain, repeat.getStart(), repeat.getEnd());
                            builder.append(pdbId + "\t" + repeat.getStart() + "-" + repeat.getEnd() + "\t" + NumberFormatUtils.format(repeat.getAvgLength()) + "\t" + repeat.getRepeats().size() + "\t" + scopeDomain + "\n");

                        }
                    }
//                    else {
//                        String scopeDomain = ScopeClassifierML.getInstance().classifyPdb(pdbCode, pdbChain, 0, 0);
//                        builder.append(pdbId + "\t" + "" + "\t" + parseDomain(scopeDomain) + "\n");
//                    }
                }
            } catch (Exception ex) {
            }

        }//end


        DataIO.writeToFile(builder.toString(), "output/trsPDB.classify.txt");

    }

    public void classifyTRsBasedOnCATH() throws Exception {
        // read tapo large-scale output
        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/resultLargeScale/tapo-v1.1.2/PDBTRs_SimilarRepeatsDB_v1.1.2.out");
        StringBuilder builder = new StringBuilder();
        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = fastaTRs.size();
        for (Map.Entry<String, TaPoFastaFormat> entry : fastaTRs.entrySet()) {
//            System.out.println(entry.getKey());
            //System.out.println(">>>>>>" + entry.getValue());
            bar.update(process, sizeOfProcess);
            process++;
            try {
                String pdbId = entry.getKey();
                String pdbCode = pdbId.substring(0, 4);
                String pdbChain = pdbId.substring(5, 6);
                TaPoFastaFormat taPoFastaFormat = entry.getValue();

                RepeatFinder repeatFinder = new RepeatFinder(pdbCode, pdbChain);
                Atom[] atoms = repeatFinder.getAtoms();

                if (taPoFastaFormat.is3DRepeat()) {
                    if (taPoFastaFormat.getRepeats().size() > 0) {

                        List<Repeat> lstRepeat = new ArrayList<Repeat>();
                        for (Repeat repeat : taPoFastaFormat.getRepeats()) {
                            if (repeat.getCluster().contains("selected")) {
                                lstRepeat.add(repeat);
                            }
                        }
                        for (Repeat repeat : lstRepeat) {

                            int seqStart = PdbTools.getResSeq(atoms, repeat.getStart());
                            int seqEnd = PdbTools.getResSeq(atoms, repeat.getEnd());
                            String region = seqStart + "-" + seqEnd;

                            ProteinCATH domain = CATHClassifier.getInstance().classifyPdb(pdbCode, pdbChain, repeat.getStart(), repeat.getEnd());
                            builder.append(pdbId + "." + repeat.getStart() + "_" + repeat.getEnd() + "\t" + region + "\t" + NumberFormatUtils.format(repeat.getAvgLength()) + "\t" + repeat.getRepeats().size() + "\t" + domain + "\n");

                        }
//                        List<CathDomain> assignedClasses = CATHClassifier.getInstance().classifyPdb(pdbCode, pdbChain, repeat.getStart(), repeat.getEnd());

                    }
                }
            } catch (Exception ex) {
            }

        }//end


        DataIO.writeToFile(builder.toString(), "output/KnownTRs.cath.txt");

    }

    public String parseDomain(String domain) {
        if (domain.length() > 1 && !domain.equals("unk")) {
            String[] scops = domain.split("\\.");
//            System.out.println(domain);
            return scops[0] + "\t" + scops[1] + "\t" + scops[2] + "\t" + scops[3];
        } else {
            return domain + "\tx\tx\tx";
        }

    }
}
