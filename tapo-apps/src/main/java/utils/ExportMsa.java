package utils;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.multalign.MutilAlign;
import org.cnrs.crbm.lib.trsfinder.Repeat;
import org.cnrs.crbm.lib.utils.ProgressBar;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Export MSA from TAPO output
 * Created by pdoviet on 10/23/2015.
 */
public class ExportMsa {


    public static void main(String[] args) throws Exception {
        new ExportMsa().exportMSA("data/phylo/tree.test.in", "data/phylo/tree.msa");

    }


    public void exportMSA(String fileIn, String fileOut) throws Exception {
        StringBuilder builder = new StringBuilder();
        List<String> pdbIds = DataIO.readLines(fileIn);
        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/resultLargeScale/tapo-v1.1.2/tapo_PDBTRs_June2015_v1.1.2.out");
        ProgressBar bar = new ProgressBar();
        int process = 1;
        int sizeOfProcess = fastaTRs.size();
        for (Map.Entry<String, TaPoFastaFormat> entry : fastaTRs.entrySet()) {
            TaPoFastaFormat taPoFastaFormat = entry.getValue();
            bar.update(process, sizeOfProcess);
            process++;
            if (taPoFastaFormat.is3DRepeat() && taPoFastaFormat.getRepeats().size() > 0) {
                List<Repeat> lstRepeat = new ArrayList<Repeat>();
                for (Repeat repeat : taPoFastaFormat.getRepeats()) {
                    if (repeat.getCluster().contains("selected")) {
                        lstRepeat.add(repeat);
                    }
                }
                for (Repeat repeat : lstRepeat) {
                    String pdbId = taPoFastaFormat.getPdbID() + "." + repeat.getStart() + "_" + repeat.getEnd();
                    if (pdbIds.contains(pdbId)) {
                        // export to MSA
                        String pdbCode = pdbId.substring(0, 4);
                        String pdbChain = pdbId.substring(5, 6);
//                        System.out.println(taPoFastaFormat.getPdbID() + "." + repeat.getStart() + "_" + repeat.getEnd());
//                        System.out.println(repeat.getRepeatRegionString());
                        MutilAlign mutilAlign = new MutilAlign(pdbCode, pdbChain, repeat.getRepeatRegionString(), true);
                        builder.append(">" + pdbId + "\n");
                        builder.append(mutilAlign.getOutput() + "\n");
                    }


                }
            }


        }


        DataIO.writeToFile(builder.toString(), fileOut);
    }


}
