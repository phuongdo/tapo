package org.cnrs.crbm.nextgen;

import org.cnrs.crbm.lib.io.DataIO;
import org.cnrs.crbm.lib.io.ReadFasta;
import org.cnrs.crbm.lib.io.TaPoFastaFormat;
import org.cnrs.crbm.lib.trsfinder.Repeat;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Created by pdoviet on 8/31/2015.
 */
public class ProcessPDB {

    public static void main(String[] args) throws Exception {


        new ProcessPDB().exportFileLargeScale();
//        new ProcessPDB().exportFastaFile();
    }


    public void exportFastaFile() {



    }


    public void exportFileLargeScale() {
        StringBuffer buffer = new StringBuffer();
        // List<String> pdbs = DataIO.readLines("data/resultLargeScale/tapo-v1.1.2/TRsRegion40.out");
        List<String> pdbs = DataIO.readLines("input/test.in");
        for (int i = 0; i < pdbs.size(); i++) {
            for (int j = i + 1; j < pdbs.size() - 1; j++) {
                buffer.append(pdbs.get(i) + "\t" + pdbs.get(j) + "\n");
            }
        }
        DataIO.writeToFile(buffer.toString(), "output/pdbpairs.in");

    }

    public void test() throws Exception {

        ReadFasta readFasta = new ReadFasta();
        Map<String, TaPoFastaFormat> fastaTRs = readFasta.readTaPoFastaFormat("data/resultLargeScale/tapo-v1.1.2/PDBTRs_NEW_v1.1.2.out");

        for (Map.Entry<String, TaPoFastaFormat> entry : fastaTRs.entrySet()) {
            TaPoFastaFormat taPoFastaFormat = entry.getValue();
            if (taPoFastaFormat.is3DRepeat() && taPoFastaFormat.getRepeats().size() > 0) {
                List<Repeat> lstRepeat = new ArrayList<Repeat>();
                for (Repeat repeat : taPoFastaFormat.getRepeats()) {
                    if (repeat.getCluster().contains("selected")) {
                        lstRepeat.add(repeat);
                    }
                }
                for (Repeat repeat : lstRepeat) {

                    System.out.println(taPoFastaFormat.getPdbID() + "." + repeat.getStart() + "_" + repeat.getEnd());

                }
            }


        }


    }


}
